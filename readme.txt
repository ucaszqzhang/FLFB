# How to run the FLFB


Step 1: Make sure you have Python3 in your computer.
## How to check and/or install Python3 version - https://realpython.com/installing-python/
## Official Python website - https://www.python.org/


Step2 : Set up a virtual enviroment with python=3.9.7
>>conda create -n FLFB python=3.9.7


Step 3: Activate the virtual environment
>>conda activate FLFB


Step 4: Download the libraries required to run FLFB model in your machine
>>pip3 install -r requirements.txt


step 5: Run FLFB on the example fasta file.
>> python3  FLFB.py  -i test.fasta  -o output_file
And you can see the result in output_file.csv 


Step 6 To run FLFB on your own fasta files
## Replace "test.fasta" with your fasta file in Step 5. You can only run one file at a time.


Step 7: After running, exit or delete your virtual environment.
exit virtual environment
>>deactivate
delete virtual environment
>> conda env remove -n FLFB
