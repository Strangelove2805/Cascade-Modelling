import numpy as np
import math
import csv
import sys


# Function for converting and reading the input file
def read_input(file, part):

    # Open the input file and check for delimiters
    with open(file, 'r') as file_contents:
        if part == 1:
            d1, d2 = "<CIRCUIT>", "</CIRCUIT>"
        elif part == 2:
            d1, d2 = "<TERMS>", "</TERMS>"
        else:
            d1, d2 = "<OUTPUT>", "</OUTPUT>"
        results = []
        for line in file_contents:
            if d1 in line:
                results = []
            elif d2 in line:
                return results
            else:
                # Check for comment markers and avoid appending the line if present
                try:
                    if (line.split())[0] == '#' or list((line.split())[0])[0] == '#':
                        continue
                    results.append(line.rstrip('\n'))
                except:
                    results.append(line.rstrip('\n'))


# Function for splitting and storing impedance information
def process_circuit(data):

    # Create a storage array with 4 columns
    storage = np.zeros((len(data), 4))

    for i in range(0, len(data)):

        # Split each node into 3 pieces of information and store accordingly
        y = data[i].split()
        one = y[0].split('=')
        two = y[1].split('=')
        thr = y[2].split('=')

        # If any of the inputs are not numbers, throw an error
        if one[1].isalpha() or two[1].isalpha() or thr[1].isalpha():
            error_output(outFile)
        storage[i][0] = int(one[1])
        storage[i][1] = int(two[1])

        # Set up the designators for component type
        switch = {
            "R": 1,
            "C": 2,
            "L": 3,
            "G": 4
        }
        storage[i][3] = switch.get(thr[0])

        # If a resistor is given in Siemens, convert to Ohms
        if storage[i][3] == 4:
            storage[i][2] = 1/float(thr[1])
            storage[i][3] = 1
        else:
            storage[i][2] = float(thr[1])

    # Sorts the node array with shunt resistors before series resistors for each node pair
    storage = storage[np.argsort(storage[:, 1])]
    storage = storage[np.argsort(storage[:, 0])]

    return storage


# Function for making the list of ABCD matrices
def ABCD_matrices(node_storage, frequency):
    matrix_list = []

    # Convert Hertzian frequency to angular frequency
    w = 2 * math.pi * frequency

    for i in range(0, len(node_storage)):

        # Calculate impedance depending on the component type
        switch = {
            1: node_storage[i][2],
            2: (1/1j) * (1 / (w * node_storage[i][2])),
            3: (1j * w * node_storage[i][2])
        }

        Z = switch.get(node_storage[i][3])

        # Create a shunt matrix
        if node_storage[i][1] == 0:
            matrix_list.append(np.array([[1, 0], [1/Z, 1]]))

        # Create a series matrix
        else:
            matrix_list.append(np.array([[1, Z], [0, 1]]))

    return matrix_list


# Function to sort the source data and frequencies
def process_source(source):

    # Store the three pieces of information
    source_storage = []
    x = source[0].split()
    y = source[1].split()
    z = source[2].split()

    # If any of the information is non-numerical, throw an error
    if (x[0].split('='))[1].isalpha() or (y[0].split('='))[1].isalpha() or (z[0].split('='))[1].isalpha():
        error_output(outFile)

    RL = int((y[0].split('='))[1])

    # Store frequency information
    is_logspaced = 0
    freq_start = float((z[0].split('='))[1])
    if list((z[0].split('='))[0])[0] == 'L':
        is_logspaced = 1

    freq_end = float((z[1].split('='))[1])
    freq_no = int((z[2].split('='))[1])

    # Create the frequency array with linear spacing
    freq = np.linspace(freq_start, freq_end, freq_no)

    # If necessary, replace the linear spaced array with a logarithmically spaced one
    if is_logspaced == 1:
        freq = freq_log(freq_start, freq_end, freq_no)

    one = x[0].split('=')
    two = x[1].split('=')
    source_storage.append(float(one[1]))

    # If source resistance is given in Siemens, convert to Ohms
    if two[0] == 'GS':
        source_storage.append(1 / float(two[1]))
    else:
        source_storage.append(float(two[1]))

    return source_storage, RL, freq


# Function to store the data related to output units
def process_output(data):

    # Create arrays for all of the necessary outputs
    values = []
    units = []
    exponents = []

    # For all of the output values, append information to the arrays
    for i in range(0, len(data)):
        y = data[i].split()
        values.append(y[0])
        try:
            units.append(y[1])
        except:

            # If the unit isn't given, fill it in
            switch = {
                "Vin":  "V",
                "Vout": "V",
                "Iin":  "A",
                "Iout": "A",
                "Pin":  "W",
                "Zout": "Ohms",
                "Pout": "W",
                "Zin":  "Ohms",
                "Av":   "L",
                "Ai":   "L"
            }

            units.append(switch.get(values[i]))

    # Check if the input values have a prefix, e.g. MegaOhms
    for i in range(0, len(units)):

        x = list(units[i])
        switch = {
            "p": "p",
            "n": "n",
            "u": "u",
            "m": "m",
            "k": "k",
            "M": "M",
            "G": "G",
            "d": "dB"
        }
        exponents.append(switch.get(x[0]))
        exponents.append(switch.get(x[0]))

    return values, units, exponents


# Function to calculate I/O impedance information
def impedances(grid, RL, ZS):

    # Store An, Bn, Cn and Dn from the product matrix
    A = grid[0][0]
    B = grid[0][1]
    C = grid[1][0]
    D = grid[1][1]

    RL = float(RL)
    Zin = (A*RL+B)/(C*RL+D)
    Zout = ((D * ZS) + B) / ((C * ZS) + A)
    Av = 1 / (A + (B * (1 / RL)))
    Ai = 1 / ((C * RL) + D)

    return Zin, Zout, Av, Ai


# Function to calculate all other circuit parameters
def parameters(VS, RL, ZS, Zin, Av, Zout):

    Iin = VS / (ZS + Zin)
    Vin = Zin * Iin
    Vout = Av * Vin
    Iout = Vout / RL
    Pin = Vin * np.conj(Iin)
    Pout = Vout * np.conj(Iout)

    return Vin, Iin, Vout, Iout, Pin, Pout


# Overarching function to create the output information grid
def calculate_answers(RL, VS, ZS, freqs, values, units, exponents):

    AnswerMatrix = []
    value_list = ["Freq"]
    unit_list = ["Hz"]

    # For all of the frequencies in the frequency list
    for i in range(0, (len(freqs) + 2)):

        # For the first loop iteration, create the value headers
        if i == 0:
            for j in range(0, len(values)):
                if exponents[(j+1)*2-2] == "dB":
                    value_list.append("|" + values[j] + "|")
                    value_list.append("/_" + values[j])
                else:
                    value_list.append("Re(" + values[j] + ")")
                    value_list.append("Im(" + values[j] + ")")
            AnswerMatrix.append(value_list)
            continue

        # For the second loop iteration, create the unit headers
        if i == 1:
            for j in range(0, len(units)):
                unit_list.append(units[j])
                unit_list.append(units[j])
            AnswerMatrix.append(unit_list)
            continue

        F = freqs[i-2]

        # Calculate the ABCD matrix and all parameters for the current frequency
        matrix_list = ABCD_matrices(node_storage, F)

        matrix_product = np.linalg.multi_dot(matrix_list)

        Zin, Zout, Av, Ai = impedances(matrix_product, RL, ZS)

        Vin, Iin, Vout, Iout, Pin, Pout = parameters(VS, RL, ZS, Zin, Av, Zout)

        parameter_container = ['{:12.3e}'.format(F)]
        for j in range(0, len(value_list)-1):

            # For a particular column header, find the output that corresponds to it
            switch = {
                "Re(Vin)":  np.real(Vin),   "Im(Vin)":  np.imag(Vin),
                "Re(Vout)": np.real(Vout),  "Im(Vout)": np.imag(Vout),
                "Re(Iin)":  np.real(Iin),   "Im(Iin)":  np.imag(Iin),
                "Re(Iout)": np.real(Iout),  "Im(Iout)": np.imag(Iout),
                "Re(Pin)":  np.real(Pin),   "Im(Pin)":  np.imag(Pin),
                "Re(Zout)": np.real(Zout),  "Im(Zout)": np.imag(Zout),
                "Re(Pout)": np.real(Pout),  "Im(Pout)": np.imag(Pout),
                "Re(Zin)":  np.real(Zin),   "Im(Zin)":  np.imag(Zin),
                "Re(Av)":   np.real(Av),    "Im(Av)":   np.imag(Av),
                "Re(Ai)":   np.real(Ai),    "Im(Ai)":   np.imag(Ai),

                # Decibel Outputs
                "|Vin|":    20*math.log(math.sqrt((np.real(Vin) ** 2) + (np.imag(Vin) ** 2))),
                "|Vout|":   20*math.log(math.sqrt((np.real(Vout) ** 2) + (np.imag(Vout) ** 2))),
                "|Iin|":    20*math.log(math.sqrt((np.real(Iin) ** 2) + (np.imag(Iin) ** 2))),
                "|Iout|":   20*math.log(math.sqrt((np.real(Iout) ** 2) + (np.imag(Iout) ** 2))),
                "|Pin|":    20*math.log(math.sqrt((np.real(Pin) ** 2) + (np.imag(Pin) ** 2))),
                "|Zout|":   20*math.log(math.sqrt((np.real(Zout) ** 2) + (np.imag(Zout) ** 2))),
                "|Pout|":   20*math.log(math.sqrt((np.real(Pout) ** 2) + (np.imag(Pout) ** 2))),
                "|Zin|":    20*math.log(math.sqrt((np.real(Zin) ** 2) + (np.imag(Zin) ** 2))),
                "|Av|":     20*math.log(math.sqrt((np.real(Av) ** 2) + (np.imag(Av) ** 2))),
                "|Ai|":     20*math.log(math.sqrt((np.real(Ai) ** 2) + (np.imag(Ai) ** 2))),
                "/_Vin":    math.atan(np.imag(Vin)/np.real(Vin)),
                "/_Vout":   math.atan(np.imag(Vout)/np.real(Vout)),
                "/_Iin":    math.atan(np.imag(Iin)/np.real(Iin)),
                "/_Iout":   math.atan(np.imag(Iout)/np.real(Iout)),
                "/_Pin":    math.atan(np.imag(Pin)/np.real(Pin)),
                "/_Zout":   math.atan(np.imag(Zout)/np.real(Zout)),
                "/_Pout":   math.atan(np.imag(Pout)/np.real(Pout)),
                "/_Zin":    math.atan(np.imag(Zin)/np.real(Zin)),
                "/_Av":     math.atan(np.imag(Av)/np.real(Av)),
                "/_Ai":     math.atan(np.imag(Ai)/np.real(Ai)),
            }

            ans = switch.get(AnswerMatrix[0][j+1])

            # If there are any prefixes for the units, multiply the out value by a correction value
            if exponents[j] is not None:
                switch = {
                    "p":    10 ** 12,
                    "n":    10 ** 9,
                    "u":    10 ** 6,
                    "m":    10 ** 3,
                    "k":    10 ** -3,
                    "M":    10 ** -6,
                    "G":    10 ** -9,
                    "dB":   1
                }

                ans *= switch.get(exponents[j])

            parameter_container.append('{:12.3e}'.format(ans))

        AnswerMatrix.append(parameter_container)

    # Correct the position of the output data with 10 characters for column 1 and 11 for the rest
    for i in range(0, len(AnswerMatrix)):
        for j in range(0, len(AnswerMatrix[0])):
            if j == 0:
                AnswerMatrix[i][j] = ((str(AnswerMatrix[i][j])).strip()).rjust(10)
            else:
                AnswerMatrix[i][j] = ((str(AnswerMatrix[i][j])).strip()).rjust(11)

    return AnswerMatrix


# Function for writing the rows of the answer matrix to a csv file
def write(AnswerMatrix, outFile):
    with open(outFile, 'w', newline='') as file:
        writer = csv.writer(file)
        for i in range(0, len(AnswerMatrix)):
            writer.writerow(AnswerMatrix[i])


# Function to produce a blank csv file if an error is detected
def error_output(outFile):
    with open(outFile, 'w') as file:
        pass
    sys.exit()


# Function to create a logarithmically spaced list of frequencies
def freq_log(Fstart, Fend, Nfreqs):
    log_start = (math.floor(math.log10(abs(Fstart))))
    log_end = (math.floor(math.log10(abs(Fend))))
    freq_log = np.logspace(log_start, log_end, Nfreqs)
    return freq_log


# Store the input / output file information from the command line
inFile = sys.argv[1]
outFile = sys.argv[2]

# Read each part of the input information
data1 = list(read_input(inFile, 1))
data2 = list(read_input(inFile, 2))
data3 = list(read_input(inFile, 3))

# Create the impedance information storage
node_storage = process_circuit(data1)

# Store source info, load resistance and frequencies
source_storage, RL, frequencies = process_source(data2)
VS = source_storage[0]
ZS = source_storage[1]

# Store output information
values, units, exponents = process_output(data3)

# Create the answer grid and write it to a csv file defined by outFile
AnswerMatrix = calculate_answers(RL, VS, ZS, frequencies, values, units, exponents)
write(AnswerMatrix, outFile)
