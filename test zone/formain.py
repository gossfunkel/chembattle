
def sayHi():
	return 'hi'
try:
	print(sayHi() for in range(10))
except SyntaxError:
	print('ya fucked it')