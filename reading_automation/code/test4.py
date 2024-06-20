import usb.core
import usb.util

# Vendor and Product IDs for the USB mouse
VENDOR_ID = 0x18F8
PRODUCT_ID = 0x1286 

# Find the USB device
device = usb.core.find(idVendor=VENDOR_ID, idProduct=PRODUCT_ID)

if device is None:
    raise ValueError("Device not found")

try:
    # Detach kernel driver if active
    if device.is_kernel_driver_active(0):
        device.detach_kernel_driver(0)

    # Set configuration
    device.set_configuration()

    # Claim interface
    usb.util.claim_interface(device, 0)

    # Read data from the mouse
    while True:
        # Read data from the endpoint (assuming endpoint 0)
        data = device.read(0x81, 8, timeout=100)
        
        # Print the raw data
        print("Raw data:", data)
        
        # Add your processing code here
        
finally:
    # Release the interface
    usb.util.release_interface(device, 0)

    # Reattach the kernel driver
    device.attach_kernel_driver(0)
