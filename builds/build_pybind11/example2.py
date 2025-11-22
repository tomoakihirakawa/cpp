#!/usr/bin/env python3
import logging
from ina226 import INA226
from time import sleep

# Configure logging to print INFO messages
logging.basicConfig(level=logging.INFO)

if __name__ == "__main__":
    try:
        # Initialize the INA226 using the default I2C address (0x40)
        # Set the maximum expected current to 25A and shunt resistance to 0.002 Ohms
        print("===================================================Begin to read")
        ina = INA226(shunt_ohms=0.002, max_expected_amps=25)
        
        print("===================================================Begin to configure")
        # Configure the INA226 with default settings
        ina.configure()
        
        print("===================================================Begin to set low battery")
        # Wake up the INA226 from sleep mode
        ina.wake()

        print("===================================================Begin to read")
        # Continuously read and print measurements
        while True:
                        
            if ina.is_conversion_ready():  # Check if conversion is ready
                print(f"Bus Voltage    : {ina.voltage():.3f} V")
                print(f"Bus Current    : {ina.current():.3f} mA")
                print(f"Supply Voltage : {ina.supply_voltage():.3f} V")
                print(f"Shunt voltage  : {ina.shunt_voltage():.3f} mV")
                print(f"Power          : {ina.power():.3f} mW")
                sleep(1)  # Delay between readings

    except (KeyboardInterrupt, Exception) as e:
        # Handle any exceptions or interrupts
        logging.error("An error occurred: %s", e)
        ina.sleep()  # Put the INA226 into sleep mode on exception
