// some functions to convert data into TeX voor MathJS

function vector_to_TeX(vector):
{
	result = "\\begin{pmatrix}"
	for element in vector[0:len(vector)-1]:
		result += str(element) + " \\\\"
	result += str(vector[len(vector)]) + " \\end{pmatrix}}

	return result
}


