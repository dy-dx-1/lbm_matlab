import pandas as pd

"""
Python program that opens an excel file in my computer with path: "D:\École\PI3\verfication_cavite.xlsx". It then takes the data from two columns 
starting at A3 and B3 respectively, but of variable length vertically, and copies relevant data into columns starting at M3 and N3 and ending at M19 and N19. 
The data relevant to copy is a map (basically a linear interpolation) of the original columns (of size Nx1 each) to the output columns at M3 and N3, 
each of size (17x1). Consider that the output columns are themselves selected data from 129 to 1, where the 17 cells represent the following indexes
starting at M3 and ending at M19: [129, 126, 125, 124, 123, 110, 96, 80, 65, 59, 37, 23, 14, 10, 9, 8, 1]. The program then saves the excel file with the new data.
"""
# Open the Excel file and specify the sheet name
file_path = r"C:\Users\Nicolas\Desktop\PI3_calcs\verification_cavite.xlsx"
sheet_name = "test"  # Replace with the actual sheet name
df = pd.read_excel(file_path, sheet_name=sheet_name)

# Get the data from columns A and B starting at A3 and B3
data_A = df.iloc[1:, 0].values
data_B = df.iloc[1:, 1].values

if len(data_A) != len(data_B): raise ValueError("The data in columns A and B must have the same length")
input_size = len(data_A) # size of the input data (num of nodes in the mesh)

# Create a map of the original columns to the output columns
output_indexes = [129, 126, 125, 124, 123, 110, 96, 80, 65, 59, 37, 23, 14, 10, 9, 8, 1]

# Create the associated output map for the indexes to exctract from the original columns, using a Lagrange interpolation 
# mapping 1 to input_size while considering lagrange functions are max and nul at 1 and 129 
output_map = [round((1*((N-129)/(1-129)))+(input_size*((N-1)/(129-1)))) for N in output_indexes]
# now rewrite output map so the indexes match the indexation, in other words 
# index N is going to be 1 and 1 is going to be N 
output_map_r = [input_size - i + 1 for i in output_map]
# Copy the relevant data to columns M and N starting at M3 and N3
df.iloc[1:18, 11] = [output_map_i for output_map_i in output_map]
df.iloc[1:18, 12] = [data_A[i-1] for i in output_map_r]
df.iloc[1:18, 13] = [data_B[i-1] for i in output_map_r]

# Now calculate the relative error between the reference data and the mapped nodes 
# So calculate error between G3 and H3 and M3 and N3, we store the result in G23 and H23 
# at cells where ref value is 0, we dont calculate and just note 0 
df.iloc[22:38, 6] = [round(100*abs((df.iloc[i, 6] - df.iloc[i, 12])/df.iloc[i, 6])) for i in range(1, 17)] # skipping last value 0
df.iloc[23:38, 7] = [round(100*abs((df.iloc[i, 7] - df.iloc[i, 13])/df.iloc[i, 7])) for i in range(2, 17)] # skipping first and last value 0 
df.iloc[38, 6] = 0
df.iloc[22, 7] = 0
df.iloc[38, 7] = 0
# Save the modified Excel file
with pd.ExcelWriter(file_path, engine='openpyxl', mode='a') as writer:
    df.to_excel(writer, sheet_name=f"{sheet_name}_", index=False)
    # delete the old sheet
    writer.book.remove(writer.book[sheet_name]) 
    # rename the new sheet to the old sheet name
    writer.book[f"{sheet_name}_"].title = sheet_name

print("Process successful") 
print("Number of nodes in the mesh: ", input_size)