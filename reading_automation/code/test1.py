from kivy.app import App
from kivy.uix.label import Label
from kivy.clock import Clock
from plyer import gyroscope
from plyer.utils import platform

# Function to request permissions
def request_permissions():
    if platform == 'android':
        from android.permissions import request_permissions, Permission
        request_permissions([Permission.BODY_SENSORS])

class GyroscopeApp(App):
    def build(self):
        self.label = Label(text="Gyroscope data will appear here")
        self.start_gyroscope()
        return self.label

    def start_gyroscope(self):
        try:
            gyroscope.enable()
            Clock.schedule_interval(self.update_gyroscope, 1.0 / 20.0)  # Update 20 times per second
        except NotImplementedError:
            self.label.text = "Gyroscope is not supported on this device."

    def update_gyroscope(self, dt):
        val = gyroscope.rotation
        if val is not None:
            self.label.text = f"Gyroscope data:\nX: {val[0]:.2f}\nY: {val[1]:.2f}\nZ: {val[2]:.2f}"
        else:
            self.label.text = "No gyroscope data available."

    def on_stop(self):
        gyroscope.disable()

if __name__ == '__main__':
    request_permissions()
    GyroscopeApp().run()
