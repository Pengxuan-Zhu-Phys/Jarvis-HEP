#!/usr/bin/env python3 
import numpy as np

# Assuming data is your numpy array
data = np.array([['0.8771738677376256', '0.09608218935151425', 'ca2f9fbe-5d91-415a-b290-d3a6dae00681'],
                 ['0.772848475187952', '0.060475599899772026', '510b4b90-710c-4e37-8fa2-96080b1ade24'],
                 ['0.9386403998656188', '0.36382096193783264', '76000ac4-a804-4d50-9afe-2d5f5818230c'],
                 ['0.278594222492493', '0.3210888643998009', '73bee864-f6d8-47e5-99dc-0933345f308f'],
                 ['0.705034293205624', '0.8655975247450612', 'd84da85a-5608-482f-8ba3-49153d05abfd'],
                 ['0.669473984645294', '0.14072049742181392', 'c565a4d2-0f01-4361-9fa1-7838ea22a2b9']])

# Split into two arrays
float_data = data[:, :-1].astype(float)  # Convert the first two columns to float
str_data = data[:, -1]  # Keep the last column as string

print(float_data, str_data)
