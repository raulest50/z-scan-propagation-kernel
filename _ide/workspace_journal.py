# 2025-07-13T16:00:53.358111700
import vitis

client = vitis.create_client()
client.set_workspace(path="z_scan_acceleration_ovr")

comp = client.get_component(name="propagation_kernel_v1")
comp.run(operation="C_SIMULATION")

comp.run(operation="C_SIMULATION")

comp.run(operation="C_SIMULATION")

comp.run(operation="C_SIMULATION")

comp.run(operation="C_SIMULATION")

comp.run(operation="C_SIMULATION")

comp.run(operation="C_SIMULATION")

comp.run(operation="C_SIMULATION")

comp.run(operation="C_SIMULATION")

comp.run(operation="C_SIMULATION")

