import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

# read to files and store them: U_matrix_1.2_pim.txt and U_matrix_1.2_us.txt

# read the files
U_matrix_1_2_pim = np.loadtxt('U_matrix_1.6_pim.txt')
U_matrix_1_2_us = np.loadtxt('U_matrix_1.6_us.txt')

absolute_error = (U_matrix_1_2_pim - U_matrix_1_2_us)
sum_error = np.sum(absolute_error)
# save the percent error to a file
np.savetxt('percent_error.txt', absolute_error, fmt='%1.4e')

print('sum_error =',sum_error)