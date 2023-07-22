# Compressibility Factor Calculation for Natural Gas using AGA8 Method
# Overview
This Python script calculates the compressibility factor of natural gas using the AGA8 method. The compressibility factor (Z-factor) is an essential parameter in the equation of state for natural gases, which allows accurate estimation of the gas behavior at different pressure and temperature conditions. The AGA8 method is widely used in the industry to determine the Z-factor of natural gas.

The script utilizes the AGA8-DC92 correlation to calculate the compressibility factor based on the gas composition, pressure, and temperature inputs. It is assumed that you have prior knowledge of the AGA8 method and are familiar with the required input parameters.
 
# Requirements
math
Ensure that you have Python installed on your system and install the necessary dependencies using the following command:
pip install math
# Usage
# 1.Input Data:

The AGA8 method requires the following input parameters:
Gas composition (mole fractions of individual components)
Pressure (in Pascals or any consistent unit)
Temperature (in Kelvin)
Make sure to have the correct gas composition data, pressure, and temperature values before running the script.

# 2.Running the Script:

Directly input it within the script. Then, run the script using the following command:

python compressibility_factor.py

# Output:
The script will calculate the compressibility factor (Z-factor) for the given gas composition, pressure, and temperature inputs. The output will be displayed in the console.

# Limitations

The script assumes that the input data is valid and correctly represents the gas composition, pressure, and temperature.
The AGA8-DC92 correlation may not be suitable for all gas mixtures. Ensure that the correlation is applicable to your specific gas composition.
The code may not cover all the aspects of error handling or input validation for production-grade applications. Enhancements can be made based on specific use cases and requirements.

# References

AGA Report No. 8 (1985), "Compressibility Factors of Natural Gas and Other Related Hydrocarbon Gases," American Gas Association.
Please feel free to modify and improve the script to suit your specific needs. If you encounter any issues or have suggestions for enhancements, feel free to contribute to the project.
