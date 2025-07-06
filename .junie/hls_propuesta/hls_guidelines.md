# Guía para Agentes AI: Generación de C++ HLS para Vitis

Esta guía está diseñada para expertos en **Xilinx Vitis HLS** con conocimientos profundos de las guías UG902 y UG1399. El objetivo es convertir especificaciones en Python en **C++ HLS** optimizado para Vitis.

## 1. Contexto y comentarios
- Antes de cada función, comenta brevemente su propósito, la latencia estimada y el uso aproximado de recursos (LUT, FF, BRAM, DSP).
- Mantén funciones pequeñas y modulares. Evita bucles anidados profundos sin pragmas.

## 2. Directivas de síntesis
Aplica siempre las siguientes pragmas para controlar paralelismo, latencia y uso de recursos:
```cpp
#pragma HLS pipeline II=<valor>
#pragma HLS unroll factor=<n>
#pragma HLS dataflow
#pragma HLS partition variable=<array> complete
```

## 3. Interfaces AXI

Define puertos de control y memoria claramente:

```cpp
#pragma HLS INTERFACE s_axilite port=return bundle=CTRL
#pragma HLS INTERFACE s_axilite port=<param> bundle=CTRL
#pragma HLS INTERFACE m_axi     port=<mem>   offset=slave bundle=G_MEM
#pragma HLS INTERFACE axis      port=<stream_in>
#pragma HLS INTERFACE axis      port=<stream_out>
```

* **AXI4-Lite** (`s_axilite`) para registros de control.
* **AXI4-Master** (`m_axi`) para accesos a DDR (burst transfers).
* **AXI-Stream** (`axis`) para flujos de datos de alta velocidad.

### 3.1 Interfaces avanzadas

Para protocolos de handshaking personalizados y control de flujo:

```cpp
// Interfaz ap_hs (handshake) para protocolos personalizados
#pragma HLS INTERFACE ap_hs port=data_in
#pragma HLS INTERFACE ap_hs port=data_out

// Interfaz ap_fifo para comunicación FIFO entre módulos
#pragma HLS INTERFACE ap_fifo port=stream_in
#pragma HLS INTERFACE ap_fifo port=stream_out

// Interfaz ap_ctrl_chain para encadenar kernels con control en serie
#pragma HLS INTERFACE ap_ctrl_chain port=return
```

Ejemplo completo de interfaz ap_ctrl_chain:

```cpp
void kernel_chain(int *in, int *out, int size) {
    #pragma HLS INTERFACE ap_ctrl_chain port=return
    #pragma HLS INTERFACE m_axi port=in bundle=gmem0
    #pragma HLS INTERFACE m_axi port=out bundle=gmem1

    // Esta interfaz permite iniciar el siguiente kernel antes de que
    // este termine completamente, mejorando el rendimiento en pipeline
    for (int i = 0; i < size; i++) {
        #pragma HLS PIPELINE II=1
        out[i] = in[i] * 2;
    }
}
```

## 4. Tipos de datos de precisión fija

Minimiza área y maximiza rendimiento ajustando el ancho de bits:

```cpp
#include <ap_int.h>
#include <ap_fixed.h>

ap_int<16>   a_int;
ap_fixed<16,6> b_fixed;
```

* `ap_int<WIDTH>` / `ap_uint<WIDTH>` para enteros de WIDTH bits.
* `ap_fixed<TOTAL,INT>` para control exacto de parte entera/fraccionaria.

## 5. Patrones de diseño modular y reutilizable

* Usa `hls::stream<T>` para conectar submódulos en modo **dataflow**.
* Encapsula lógica común (filtros, convoluciones) en plantillas C++.
* Emplea `ap_ctrl_chain` para encadenar kernels con control en serie.

## 6. Precisión fija y análisis de bit-width

* Realiza una fase de **profiling** con vectores de prueba para determinar bit-widths mínimos.
* Documenta cada elección de bit-width en comentarios con su impacto estimado en lógica y latencia.

## 7. Optimización de memoria y acceso a datos

* Aplica `#pragma HLS ARRAY_PARTITION` para eliminar cuellos de botella de memoria.
* Configura y verifica **burst transfers** en el log de síntesis para accesos óptimos.
* Analiza patrones de acceso y ajusta tamaños de ráfaga.

## 8. Integración de CI/CD

* Configura un pipeline (GitHub Actions, GitLab CI…) que ejecute `csim`, `csynth`, `cosim` y `export_ip` en cada PR.
* Incluye linters de estilo C++ y comprobaciones de presencia de pragmas.
* Reporta métricas clave (II, latencia, recursos) en la interfaz de revisión.

### 8.1 Ejemplo de GitHub Actions para HLS

```yaml
# .github/workflows/hls_ci.yml
name: HLS CI

on:
  push:
    branches: [ main ]
  pull_request:
    branches: [ main ]

jobs:
  hls_build:
    runs-on: self-hosted
    steps:
      - uses: actions/checkout@v3

      - name: Setup Vitis
        run: |
          source /opt/Xilinx/Vitis/2025.1/settings64.sh
          echo "XILINX_VITIS set to $XILINX_VITIS"

      - name: C++ Lint
        run: |
          cpplint --recursive --filter=-legal/copyright,-build/include_subdir src/

      - name: C Simulation
        run: |
          cd hls_project
          vitis_hls -f scripts/run_csim.tcl

      - name: C Synthesis
        run: |
          cd hls_project
          vitis_hls -f scripts/run_csynth.tcl

      - name: Extract Metrics
        run: |
          cd hls_project
          python scripts/parse_report.py --report_file syn/report/csynth.xml

      - name: Co-Simulation
        run: |
          cd hls_project
          vitis_hls -f scripts/run_cosim.tcl

      - name: Export IP
        run: |
          cd hls_project
          vitis_hls -f scripts/export_ip.tcl

      - name: Archive Results
        uses: actions/upload-artifact@v3
        with:
          name: hls-artifacts
          path: |
            hls_project/syn/report/
            hls_project/impl/ip/
```

### 8.2 Script para extraer métricas de síntesis

```python
# scripts/parse_report.py
import argparse
import xml.etree.ElementTree as ET
import json

def parse_hls_report(report_file):
    tree = ET.parse(report_file)
    root = tree.getroot()

    # Extraer métricas de latencia
    latency = root.find('.//PerformanceEstimates/SummaryOfOverallLatency')
    latency_min = latency.find('Best-caseLatency').text
    latency_max = latency.find('Worst-caseLatency').text

    # Extraer métricas de recursos
    resources = root.find('.//AreaEstimates/Resources')
    lut = resources.find('LUT').text
    ff = resources.find('FF').text
    bram = resources.find('BRAM_18K').text
    dsp = resources.find('DSP').text

    # Extraer métricas de timing
    timing = root.find('.//AreaEstimates/AvailableResources')
    target_cp = timing.find('AVAIL_CLK').text

    # Crear diccionario de resultados
    results = {
        "latency": {
            "min": latency_min,
            "max": latency_max
        },
        "resources": {
            "LUT": lut,
            "FF": ff,
            "BRAM": bram,
            "DSP": dsp
        },
        "timing": {
            "target_cp": target_cp
        }
    }

    # Extraer información de loops
    loops = root.findall('.//LoopLatency')
    loop_info = []

    for loop in loops:
        loop_name = loop.find('Name').text
        loop_ii = loop.find('InterIterationInterval').text
        loop_info.append({
            "name": loop_name,
            "II": loop_ii
        })

    results["loops"] = loop_info

    return results

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Parse HLS synthesis report')
    parser.add_argument('--report_file', type=str, required=True, help='Path to csynth.xml report file')
    parser.add_argument('--output', type=str, default='metrics.json', help='Output JSON file')

    args = parser.parse_args()

    results = parse_hls_report(args.report_file)

    # Imprimir resultados en consola
    print("=== HLS Synthesis Results ===")
    print(f"Latency (min/max): {results['latency']['min']}/{results['latency']['max']} cycles")
    print(f"Resources: LUT={results['resources']['LUT']}, FF={results['resources']['FF']}, BRAM={results['resources']['BRAM']}, DSP={results['resources']['DSP']}")

    print("\nLoop Performance:")
    for loop in results['loops']:
        print(f"  {loop['name']}: II={loop['II']}")

    # Guardar resultados en JSON
    with open(args.output, 'w') as f:
        json.dump(results, f, indent=2)

    print(f"\nResults saved to {args.output}")
```

## 9. Automatización y scripting con Python en Vitis

* Usa la **API Python** de Vitis para invocar fases de HLS en batch mode.
* Captura comandos GUI y transpórtalos a un script `.py` reproducible.
* Permite pasar parámetros dinámicos (N, unroll factor) como argumentos al script.

## 10. Directrices de pruebas y validación

* Genera un **testbench** en C++ que cubra casos límite y realice comparaciones bit a bit.
* Compara salidas de `csim` con un "golden model" en Python o MATLAB antes de `csynth`.
* Itera: tras cada síntesis, ajusta pragmas para cumplir objetivos de II y latencia, documentando el impacto.

## 11. Ejemplos de la documentación oficial AMD Xilinx

### 11.1 Pipeline de bucles (UG902)

```cpp
void foo(char x, char a, char b, char c) {
  #pragma HLS pipeline II=1
  char y = x*a + b + c;
}
```

### 11.2 Unrolling de bucles (UG902)

```cpp
for (int i = 0; i < N; i++) {
  #pragma HLS unroll factor=4
  y[i] = a * x[i] + y[i];
}
```

### 11.3 Dataflow entre funciones (UG902)

```cpp
#pragma HLS dataflow
read_stream(in, buf);
process(buf, tmp);
write_stream(tmp, out);
```

### 11.4 Interfaces AXI (UG902)

```cpp
#pragma HLS INTERFACE s_axilite port=return bundle=CTRL
#pragma HLS INTERFACE m_axi     port=in1  offset=slave bundle=G_MEM
#pragma HLS INTERFACE m_axi     port=in2  offset=slave bundle=G_MEM
#pragma HLS INTERFACE m_axi     port=out  offset=slave bundle=G_MEM
#pragma HLS INTERFACE axis      port=stream_in
#pragma HLS INTERFACE axis      port=stream_out
```

### 11.5 Tipos de datos de precisión fija (UG902)

```cpp
#include <ap_int.h>
#include <ap_fixed.h>

ap_int<16> a_int;
ap_fixed<16,6> b_fixed;
```

### 11.6 Acceso en ráfaga (UG1399)

```cpp
void vadd(const int *in1, const int *in2, int *out, int size) {
  #pragma HLS INTERFACE m_axi port=in1  offset=slave bundle=G_MEM
  #pragma HLS INTERFACE m_axi port=in2  offset=slave bundle=G_MEM
  #pragma HLS INTERFACE m_axi port=out  offset=slave bundle=G_MEM
  for (int i = 0; i < size; i++) {
    #pragma HLS pipeline II=1
    out[i] = in1[i] + in2[i];
  }
}
```

### 11.7 Modularización con `hls::stream` (UG902)

```cpp
#include <hls_stream.h>

void stage1(hls::stream<int>& in, hls::stream<int>& mid) {
  #pragma HLS pipeline II=1
  int v = in.read();
  mid.write(v * 2);
}

void stage2(hls::stream<int>& mid, hls::stream<int>& out) {
  #pragma HLS pipeline II=1
  int v = mid.read();
  out.write(v + 1);
}

void top(hls::stream<int>& in, hls::stream<int>& out) {
  hls::stream<int> mid;
  #pragma HLS dataflow
  stage1(in, mid);
  stage2(mid, out);
}
```

### 11.8 Ejemplo de "Perfect Loop" (UG902)

```cpp
void perfect(int in[3], int out[3]) {
  #pragma HLS pipeline II=1
  for (int i = 0; i < 3; i++) {
    out[i] = in[i] * 2;
  }
}
```

### 11.9 Testbench y verificación (UG902)

```cpp
#include "saxpy.h"
#include <iostream>

int main() {
  const int N = 16;
  float a = 2.5f;
  float x[N], y[N], y_ref[N];
  for (int i = 0; i < N; i++) {
    x[i] = i;
    y[i] = i;
    y_ref[i] = a * x[i] + y[i];
  }
  saxpy(a, x, y, N);
  for (int i = 0; i < N; i++) {
    if (y[i] != y_ref[i]) {
      std::cerr << "Error en i=" << i << "\n";
      return 1;
    }
  }
  std::cout << "Test exitoso\n";
  return 0;
}
```

## 12. Equivalencias Python a C++ HLS

### 12.1 Tipos de datos

| Python | C++ HLS | Notas |
|--------|---------|-------|
| `int` | `int`, `ap_int<W>` | Usar `ap_int<W>` para optimizar ancho de bits |
| `float` | `float`, `ap_fixed<W,I>` | Usar `ap_fixed<W,I>` para precisión fija |
| `bool` | `bool` | Mismo comportamiento |
| `list`, `array` | Arrays estáticos, `hls::vector<T,N>` | Tamaño debe ser conocido en tiempo de compilación |
| `numpy.ndarray` | Arrays multidimensionales | Declarar como `T array[DIM1][DIM2]` |
| `dict` | No soportado directamente | Implementar como arrays paralelos o estructuras |

### 12.2 Estructuras de control

| Python | C++ HLS | Consideraciones HLS |
|--------|---------|---------------------|
| `for i in range(n)` | `for (int i = 0; i < n; i++)` | Añadir `#pragma HLS pipeline II=1` para paralelismo |
| `while condition` | `while (condition)` | Evitar condiciones de salida impredecibles |
| `if-elif-else` | `if-else if-else` | Minimizar ramas condicionales complejas |
| List comprehension | Bucles explícitos | Desenrollar con `#pragma HLS unroll` |
| Funciones lambda | Functores o funciones inline | Preferir funciones explícitas para mejor síntesis |

## 13. Patrones de conversión Python a C++ HLS

### 13.1 Procesamiento de arrays NumPy

**Python:**
```python
def vector_add(a, b):
    return a + b  # Operación vectorizada
```

**C++ HLS:**
```cpp
template<int N>
void vector_add(const float a[N], const float b[N], float result[N]) {
    #pragma HLS INTERFACE m_axi port=a bundle=gmem0
    #pragma HLS INTERFACE m_axi port=b bundle=gmem1
    #pragma HLS INTERFACE m_axi port=result bundle=gmem2

    vector_add_kernel:
    for (int i = 0; i < N; i++) {
        #pragma HLS PIPELINE II=1
        result[i] = a[i] + b[i];
    }
}
```

### 13.2 Convolución 2D

**Python:**
```python
def conv2d(input_data, kernel):
    result = np.zeros_like(input_data)
    k_height, k_width = kernel.shape
    for i in range(input_data.shape[0] - k_height + 1):
        for j in range(input_data.shape[1] - k_width + 1):
            result[i, j] = np.sum(input_data[i:i+k_height, j:j+k_width] * kernel)
    return result
```

**C++ HLS:**
```cpp
template<int HEIGHT, int WIDTH, int K_HEIGHT, int K_WIDTH>
void conv2d(const float input[HEIGHT][WIDTH], 
            const float kernel[K_HEIGHT][K_WIDTH],
            float output[HEIGHT-K_HEIGHT+1][WIDTH-K_WIDTH+1]) {
    #pragma HLS INTERFACE m_axi port=input offset=slave bundle=gmem1
    #pragma HLS INTERFACE m_axi port=kernel offset=slave bundle=gmem2
    #pragma HLS INTERFACE m_axi port=output offset=slave bundle=gmem3

    // Buffers locales para optimizar acceso a memoria
    float kernel_local[K_HEIGHT][K_WIDTH];
    #pragma HLS ARRAY_PARTITION variable=kernel_local complete dim=0

    // Cargar kernel en memoria local
    kernel_load: for (int i = 0; i < K_HEIGHT; i++) {
        for (int j = 0; j < K_WIDTH; j++) {
            #pragma HLS PIPELINE II=1
            kernel_local[i][j] = kernel[i][j];
        }
    }

    // Ventana deslizante para convolución
    row_loop: for (int i = 0; i < HEIGHT-K_HEIGHT+1; i++) {
        col_loop: for (int j = 0; j < WIDTH-K_WIDTH+1; j++) {
            #pragma HLS PIPELINE II=1

            float sum = 0.0f;
            k_row: for (int ki = 0; ki < K_HEIGHT; ki++) {
                k_col: for (int kj = 0; kj < K_WIDTH; kj++) {
                    #pragma HLS UNROLL
                    sum += input[i+ki][j+kj] * kernel_local[ki][kj];
                }
            }
            output[i][j] = sum;
        }
    }
}
```

## 14. Optimizaciones específicas de Vitis 2025.1

### 14.1 Inferencia de interfaces mejorada

Vitis 2025.1 simplifica aún más la inferencia automática de interfaces:

```cpp
// Vitis 2025.1 infiere automáticamente interfaces AXI adecuadas
// con configuración optimizada basada en patrones de acceso
void top_function(int *input, int *output, int size) {
    // El compilador infiere:
    // #pragma HLS INTERFACE m_axi port=input burst=16 num_read_outstanding=4
    // #pragma HLS INTERFACE m_axi port=output burst=16 num_write_outstanding=4
    // #pragma HLS INTERFACE s_axilite port=size
    // #pragma HLS INTERFACE s_axilite port=return

    for (int i = 0; i < size; i++) {
        #pragma HLS PIPELINE II=1
        output[i] = input[i] * 2;
    }
}
```

### 14.2 Soporte mejorado para tipos de datos C++

```cpp
#include <vector>
#include <array>
#include <algorithm>

void modern_cpp_kernel(std::array<float, 1024>& input, std::array<float, 1024>& output) {
    #pragma HLS INTERFACE m_axi port=input
    #pragma HLS INTERFACE m_axi port=output

    std::transform(input.begin(), input.end(), output.begin(), 
                  [](float x) { return x * 2.0f; });
}
```

## 15. Integración con Python mediante pybind11

Para facilitar la verificación entre modelos Python y C++ HLS:

```cpp
// En archivo binding.cpp
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include "hls_kernel.h"

namespace py = pybind11;

// Wrapper para llamar al kernel HLS desde Python
py::array_t<float> hls_kernel_wrapper(py::array_t<float> input) {
    auto buf_in = input.request();
    auto result = py::array_t<float>(buf_in.size);
    auto buf_out = result.request();

    float* ptr_in = static_cast<float*>(buf_in.ptr);
    float* ptr_out = static_cast<float*>(buf_out.ptr);

    // Llamar a la función HLS
    hls_kernel(ptr_in, ptr_out, buf_in.size);

    return result;
}

PYBIND11_MODULE(hls_module, m) {
    m.def("process", &hls_kernel_wrapper, "HLS kernel function");
}
```

## 16. Depuración y perfilado

### 16.1 Instrumentación para depuración

```cpp
#ifndef __SYNTHESIS__
#include <iostream>
#define DEBUG(x) std::cout << x << std::endl
#else
#define DEBUG(x)
#endif

void kernel_function(int *data, int size) {
    for (int i = 0; i < size; i++) {
        #pragma HLS PIPELINE II=1
        data[i] *= 2;
        DEBUG("Processed element " << i << ": " << data[i]);
    }
}
```

### 16.2 Análisis de rendimiento con pragmas de instrumentación

```cpp
void compute_kernel(float *input, float *output, int size) {
    #pragma HLS INTERFACE m_axi port=input
    #pragma HLS INTERFACE m_axi port=output

    #pragma HLS PERFORMANCE target_latency=1000

    for (int i = 0; i < size; i++) {
        #pragma HLS PIPELINE II=1
        output[i] = complex_function(input[i]);
    }
}
```

## 17. Patrones para algoritmos de aprendizaje automático

### 17.1 Implementación de capa fully-connected

**Python (PyTorch):**
```python
def forward(self, x):
    return torch.matmul(x, self.weights.t()) + self.bias
```

**C++ HLS:**
```cpp
template<int IN_SIZE, int OUT_SIZE>
void fc_layer(const float input[IN_SIZE], 
              const float weights[OUT_SIZE][IN_SIZE],
              const float bias[OUT_SIZE],
              float output[OUT_SIZE]) {
    #pragma HLS ARRAY_PARTITION variable=weights cyclic factor=16 dim=2

    for (int o = 0; o < OUT_SIZE; o++) {
        #pragma HLS PIPELINE II=1
        float sum = bias[o];
        for (int i = 0; i < IN_SIZE; i++) {
            #pragma HLS UNROLL factor=16
            sum += input[i] * weights[o][i];
        }
        output[o] = sum;
    }
}
```

### 17.2 Implementación de activación ReLU

**Python:**
```python
def relu(x):
    return np.maximum(0, x)
```

**C++ HLS:**
```cpp
template<int SIZE>
void relu(const float input[SIZE], float output[SIZE]) {
    for (int i = 0; i < SIZE; i++) {
        #pragma HLS PIPELINE II=1
        output[i] = (input[i] > 0) ? input[i] : 0;
    }
}
```

## 18. Consideraciones de memoria y ancho de banda

### 18.1 Optimización de patrones de acceso a memoria

```cpp
void optimize_memory_access(const float *input, float *output, int rows, int cols) {
    float buffer[MAX_ROWS][MAX_COLS];
    #pragma HLS ARRAY_PARTITION variable=buffer cyclic factor=16 dim=2

    // Cargar datos con acceso eficiente en ráfaga
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            #pragma HLS PIPELINE II=1
            buffer[i][j] = input[i*cols + j];
        }
    }

    // Procesar datos desde buffer local
    process_data:
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            #pragma HLS PIPELINE II=1
            // Procesamiento...
            output[i*cols + j] = buffer[i][j] * 2;
        }
    }
}
```

### 18.2 Uso de memoria en chip vs. externa

```cpp
void memory_hierarchy_example(const float *input, float *output, int size) {
    // BRAM para datos frecuentemente accedidos
    float local_buffer[1024];
    #pragma HLS ARRAY_PARTITION variable=local_buffer cyclic factor=4

    // URAM para conjuntos de datos más grandes
    float large_buffer[4096];
    #pragma HLS RESOURCE variable=large_buffer core=RAM_2P_URAM

    // Cargar datos críticos en BRAM
    for (int i = 0; i < 1024; i++) {
        #pragma HLS PIPELINE II=1
        local_buffer[i] = input[i];
    }

    // Cargar datos secundarios en URAM
    for (int i = 0; i < 4096; i++) {
        #pragma HLS PIPELINE II=1
        large_buffer[i] = input[1024 + i];
    }

    // Procesamiento usando jerarquía de memoria optimizada
    // ...
}
```

---

Esta guía sirve como **plantilla única** para configurar cualquier agente de IA o entorno de codificación, garantizando que genere C++ HLS de alta calidad acorde a las mejores prácticas de AMD Xilinx.
