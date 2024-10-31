/* some functions to convert data into TeX voor MathJS
*/

function vector_to_TeX(vector)
{
	result = "\\begin{pmatrix}"
	let i = 0;
	for ( ; i<vector.length-1; i++)
		result += vector[i].toString() + " \\\\"
	result += vector[i].toString() + " \\end{pmatrix}"

	return result
}


