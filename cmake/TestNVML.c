#include <nvml.h>

int main()
{
    // This program will only link correctly if the NVML library version in the display driver
    // supports the CUDA toolkit version.
    nvmlInit();
    unsigned int nvml_device_count;
    nvmlDeviceGetCount(&nvml_device_count);
    nvmlDevice_t nvml_device_id;
    nvmlDeviceGetHandleByIndex(0, &nvml_device_id);
    nvmlPciInfo_t nvml_pci_info;
    nvmlDeviceGetPciInfo(nvml_device_id, &nvml_pci_info);
    nvmlShutdown();
    return 0;
}
