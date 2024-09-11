Code usage:
- To run code
    - python3 conservation.py
- Besides the comments describing the code...
    - read_fna_sequences is for creating a list of sequences from a txt file
    - calculate_conservation calculates rates at each position of sequences 
        - uncomment the early return statement within the function and edit the unpacking in main
        to have original conservation rates returned. Other wise leave as is and get returned a list of 
        smoothed conserved sequences and a data frame containign said smoothed sequences. 
        - I use the smoothing package in pandas and base it off of the last 90 conserved values, so the 
        first 90 postions are rendered nan in order to smooth the entirety of the data set
    - variable_regions is for finding the variable regions within the data set. they are stored as a list of tuples
        - the threshold used is the average of all the data points. Found that there were 8 variable regions
    - The other 2 output functions are for writing answers to txt files for questions 1 and 3. 