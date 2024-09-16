import usb.core
import usb.util
import time

VENDOR_ID = 0x0922
PRODUCT_ID = 0x8003
MAX_ATTEMPTS = 10


def find_device():
    device = usb.core.find(idVendor=VENDOR_ID, idProduct=PRODUCT_ID)
    while device is None:
        print("Searching for device...")
        device = usb.core.find(idVendor=VENDOR_ID, idProduct=PRODUCT_ID)
    return device

def read_data(device, endpoint):
    attempts = MAX_ATTEMPTS
    while attempts > 0:
        try:
            data = device.read(endpoint.bEndpointAddress, endpoint.wMaxPacketSize)
            grams = data[4] + (256 * data[5])
            print(str(grams) + "g")
            return
        except usb.core.USBError as e:
            device.set_configuration()
            attempts -= 1
            print("Failure! Attempts left:", attempts)

        time.sleep(1)

    print("Failed to connect")

device = find_device()
device.set_configuration()
endpoint = device[0][(0,0)][0]

while True:
    read_data(device, endpoint)
