class NumberError(Exception):
    def __init__(self):
        print("\nException: NumberError:\n the number you have typed is wrong, you have to choose between '1', '2' or '3'\n")
        
class NumberError2(Exception):
    def __init__(self):
        print("\nException: NumberError:\n the number you have typed is wrong, you have to choose between '1', '2', '3', '4' or '5'\n")
        
        
class ValueNegative(Exception):
	def __init__(self):
		print("\nException: ValueNegative:\n you have typed a negative number, in order to trim you must give a positive integer\n")
		
	
class TaskError(Exception):
	def __init__(self):
		print("\nException: the current task failed\n")
		
		
class Kmer(Exception):
	def __init__(self):
		print("\nException: the number has to be an integer between 20 and 50\n")
		
		