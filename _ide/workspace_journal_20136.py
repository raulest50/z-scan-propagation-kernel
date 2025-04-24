# 2025-04-07T16:56:16.580124700
import vitis

client = vitis.create_client()
client.set_workspace(path="z_scan_acceleration_ovr")

status = client.add_platform_repos(platform=["c:\Vws\base_platform_component"])

comp = client.create_hls_component(name = "propagation_kernel_v1",cfg_file = ["hls_config.cfg"],template = "empty_hls_component")

vitis.dispose()

