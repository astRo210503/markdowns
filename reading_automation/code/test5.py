import usb.core
import usb.util
import csv
import time

# Find the USB mouse by HID device class
device = usb.core.find( idVendor=0x18F8, idProduct=0x1286)

ep = device[0].interfaces()[0].endpoints()[0]
i =  device[0].interfaces()[0].bInterfaceNumber
device.reset()


    
device.set_configuration()
eaddr = ep.bEndpointAddress

r = device.read(eaddr,1024)
print(len(r))







# print("ep :" , ep)



