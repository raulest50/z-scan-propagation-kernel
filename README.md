# Z-Scan Acceleration Kernel

This repository contains a hardware accelerated implementation of a Z-Scan propagation kernel targeting the AMD-Xilinx KV260 board. The code was developed using the **Vitis Unified IDE 2024.2** and includes both the high level synthesis (HLS) sources and reference test data.

## Repository Layout

- `propagation_kernel_v1/` – main HLS component
  - `bpm_prop_accel.cpp` – top level kernel
  - `step_propagators.cpp`/`step_propagators.h` – propagation routines
  - `step_propagators_tb.cpp` – standalone test bench
  - `validationData/` – input and reference output files for the tests
  - `hls_config.cfg` – configuration passed to Vitis HLS
  - `vitis-comp.json` – component description used by the IDE
- `KV260_resources.md` – notes about KV260 FPGA resources

## Prerequisites

- **Vitis / Vitis HLS 2024.2** or newer installed and added to your `PATH`.
- Access to the headers shipped with Vitis (`hls_x_complex.h`, `hls_math.h`).
- This repository cloned in a writable location so the build and simulation can create intermediate files.

## Building with Vitis HLS

The kernel is built through the HLS flow using the `hls_config.cfg` file. The easiest way is via the Vitis IDE:

1. Launch *Vitis Unified IDE 2024.2*.
2. Import the component by opening `vitis-comp.json`.
3. Run **C Simulation** or **C Synthesis**. The IDE reads all options from `hls_config.cfg` and places generated results in `propagation_kernel_v1/bpm_prop_accel`.

From the command line you can perform the same step with:

```bash
vitis_hls -i vitis-comp.json
```

This command invokes the HLS tool using the configuration in `hls_config.cfg`.

## Running the Test Benches

`step_propagators_tb.cpp` provides a simple test that compares the HLS results against the reference data in `validationData/`.

To execute it inside Vitis HLS:

1. Open the component as described above.
2. Select `step_propagators_tb.cpp` as the test bench and run **C Simulation**.
3. Ensure the working directory is the repository root so the files under `propagation_kernel_v1/validationData` are found.
4. On success the simulation prints an RMS error similar to:

```
RMS error: 1.0e-03
```

## Tool Versions

The project was created with:

- **Vitis Unified IDE:** 2024.2
- **Target device:** `xck26-sfvc784-2LV-c`
- The host environment here uses `g++` version `g++ (Ubuntu 13.3.0-6ubuntu2~24.04) 13.3.0`.

