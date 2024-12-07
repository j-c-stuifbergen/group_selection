Array.prototype.indexOfMaxAbs = function(minPosition=0, upperPosition=null) {

    if (null == upperPosition) {
		upperPosition = this.length
	}

	if (this.length <= minPosition) {
        return -1;
    }
	if (this.length < upperPosition) {
        return -2;
    }
	
    var maxValue = this[minPosition];
    var maxIndex = minPosition;

    for (var i = minPosition +1; (i < this.length) && (i<upperPosition); i++) {
        if ( maxValue < Math.abs(this[i])) {
            maxIndex = i;
            maxValue = Math.abs(this[i]);
        }
    }
    return maxIndex;
}

Array.prototype.timesScalar = function(scalar)
{	var result = new Array(this.length)
	for(let i = 0; i<this.length ; i++)
	{	result[i] = this[i]*scalar
	}
	return result
}
Array.prototype.add = function(vector)
{	
	var result = new Array(this.length)
	if (this.length != vector.length)
	{	throw ("addition of vectors of different length")
	}
	for (i = 0; i< this.length; i++)
	{	result[i] = this[i] + vector[i]
	}
	return result
}
Array.prototype.subtract = function(vector)
{
	var result = new Array(this.length)
	if (this.length != vector.length)
	{	throw ("addition of vectors of different length")
	}
	for (i = 0; i< this.length; i++)
	{	result[i] = this[i] - vector[i]
	}
	return result
}
Array.prototype.inner = function(vector)
{
	var result = 0
	if (this.length != vector.length)
	{	throw ("inner product of vectors of different length")
	}
	for (i = 0; i< this.length; i++)
	{	result += this[i] * vector[i]
	}
	return result
}
Array.prototype.multiplyByMatrix = function(M)
{
	var result = Array(M.length)
	for (let i = 0; i< M.length ; i++)
	{	result[i] = M[i].inner(this)
	}
	return result
}

function findMaxAbsElement(Matrix, minRow = 0, upperRow =null, minCol = 0, upperCol = null)
{
	if (null == upperRow)
	{	upperRow = Matrix.length
	}
	if (Matrix.length < upperRow)
	{	throw ("row limit is too big for the matrix.")
	}
	// else
	if (null == upperCol)
	{
		upperCol = Matrix[0].length
	}

	let maxFound = 0
	rowIndex = minRow
	colIndex = minRow
	for ( let i = minRow; i<upperRow ; i++)
	{
		if (Matrix[i].length < upperCol)
		{	throw ("inconsistent Row lengths")
			return
		}
		highestInRow = Matrix[i].indexOfMaxAbs(minCol, upperCol)
		if ( maxFound < Math.abs(Matrix[i][highestInRow]))
		{	maxFound = Math.abs(Matrix[i][highestInRow])
		    colIndex = highestInRow
			rowIndex = i
		}
	}	
	return { rowIndex: rowIndex , columnIndex: colIndex}	
}

// Function to solve a matrix equation A * X = B
// There may be not be more columns than rows.
// The matrix can be singular, but if A X = B has no solution, an error is thrown.
// The indexes are A[row][column].
// X.length = number of columns, B.length = number of rows = A.length
// epsilon should depend on the machine precision, and perhaps on the norm of the matrix.
// When pivoting a singular matrix, elements under the last relevant row should be smaller than epsilon.
// Corresponding elements of B should be smaller than margin * epsilon.
function solveMatrixEquation(A, B, epsilon = 1e-12, margin =5)
{
  const nRows = A.length;
  var dimension = 0 // number of independent vectors
  if (0<nRows)
  {	var nCols = A[0].length
	dimension = nCols 
  }
  let colPerms = Array.from(Array(nCols).keys())
  // Augment matrix A with column vector B
  let augmentedMatrix = A.map((row, i) => [...row, B[i]]);
  
  // Forward elimination with partial pivoting
  for (let col = 0; col < nCols; col++) {
    // find the pivot element
	let pivotOrdinates = findMaxAbsElement(augmentedMatrix, col, nRows, col, nCols)
	let pivotRow = pivotOrdinates.rowIndex
	let pivotColumn = pivotOrdinates.columnIndex
    // Swap the current row with the pivot row (this doesn't change the solution)
    if (pivotRow !== col) {
      let temp = augmentedMatrix[col];
      augmentedMatrix[col] = augmentedMatrix[pivotRow];
      augmentedMatrix[pivotRow] = temp;
    }
	if (pivotColumn !== col)
	{	for (let i = 0; i<nRows; i++)
		{	temp = augmentedMatrix[i][col]
			augmentedMatrix[i][col] = augmentedMatrix[i][pivotColumn]
			augmentedMatrix[i][pivotColumn] = temp
		}
		// store the column permutations
		{  let temp = colPerms[ col] 
			colPerms[ col] = colPerms[pivotColumn]
			colPerms[pivotColumn] = temp
		}
	}

	if (epsilon < Math.abs(augmentedMatrix[col][col]))
	{
	    // Perform elimination to make all values below the pivot zero
	    for (let row = col + 1; row < nRows; row++) {
	      const factor = augmentedMatrix[row][col] / augmentedMatrix[col][col];
		  augmentedMatrix[row][col] = 0 // if not zero, it's a round-off error
	      for (let j = col+1; j <= nCols; j++) {
			augmentedMatrix[row][j] -= factor * augmentedMatrix[col][j];
	      }
	    }
	}
	else // all remaining elements are smaller than epsilon
	{
	    dimension = col // the previous row was the last independent row
	    console.log("for this matrix of size "+nRows+", there are only "+col+" independent columns")
		// all remaining elements will be ignored in back-substitution
		// So there is no more pivoting.
	    col = nCols // leave the loop
	}
  }

	// if the equation has a valid solution, all remaining elements in the last column are zero.
	// We must check, because there is no back-substitution for those elements
	for (let row = dimension ; row < nRows; row++) {
		if (Math.abs(augmentedMatrix[row][nCols]) > epsilon * margin)
		{	throw ("error: singular matrix,  no solution: After elimination, vector element ["+row+"] is "+augmentedMatrix[row][nCols]+ " > "+epsilon * margin+" =epsilon * margin")
		}
	}

  // Back substitution
  let X = new Array(nCols).fill(0);
  for (let row = dimension - 1; row >= 0; row--) {
    let sum = augmentedMatrix[row][nCols]; // Right-hand side of the equation
    for (let col = row + 1; col < dimension; col++) {
      sum -= augmentedMatrix[row][col] * X[col];
    }
    if ( 0!=sum) // then 0!=augmentedMatrix[row][row])
    X[row] = sum / augmentedMatrix[row][row];
  }
  // correct for the permutation of the columns
	var result = Array(nCols)
	for (i = 0; i<nCols ; i++)
	{  result [colPerms[i]] = X[i]
	}
  return result;
}

// find independent rows of matrix M
function findBasis(matrix, epsilon = 1e-12)
{
  const nRows = matrix.length;
  var dimension = 0 // number of independent vectors
  if (0<nRows)
  {	var nCols = A[0].length
	dimension = nCols 
  }
  let rowPerms = Array.from(Array(nRows).keys())
  
  // Forward elimination with partial pivoting
  for (let col = 0; col < nCols; col++) {
    // find the pivot element
	let pivotOrdinates = findMaxAbsElement(matrix, col, nRows, col, nCols)
	let pivotRow = pivotOrdinates.rowIndex
	let pivotColumn = pivotOrdinates.columnIndex
    // Swap the current row with the pivot row (this doesn't change the solution)
    if (pivotRow !== col) {
      let temp = matrix[col];
      matrix[col] = matrix[pivotRow];
      matrix[pivotRow] = temp;
	// store the row permutations
      temp = rowPerms[ col] 
      rowPerms[ col] = rowPerms[pivotRow]
      rowPerms[pivotRow] = temp
    }

	if (pivotColumn !== col)
	{	for (let i = 0; i<nRows; i++)
		{	temp = matrix[i][col]
			matrix[i][col] = matrix[i][pivotColumn]
			matrix[i][pivotColumn] = temp
		}
	}

	if (epsilon < Math.abs(matrix[col][col]))
	{
	    // Perform elimination to make all values below the pivot zero
	    for (let row = col + 1; row < nRows; row++) {
	      const factor = matrix[row][col] / matrix[col][col];
		  matrix[row][col] = 0 // if not zero, it's a round-off error
	      for (let j = col+1; j < nCols; j++) {
			matrix[row][j] -= factor * matrix[col][j];
	      }
	    }
	}
	else // all remaining elements are smaller than epsilon
	{
	    dimension = col // the previous row was the last independent row
	    console.log("for this matrix of size "+nRows+", there are only "+col+" independent columns")
		// all remaining elements will be ignored in back-substitution
		// So there is no more pivoting.
	    col = nCols // leave the loop
	}
  }

  // correct for the permutation of the columns
	rowPerms.splice(dimension)

	console.log("dimension is "+dimension+", the independent vectors have indexes "+rowPerms)
  return rowPerms;
}

function transpose(M)
{	// number of columns of the transpose
	nCols = M.length // = number of rows of M
	if (0 == nCols)
	{	return []
	}
	// else
	nRows = M[0].length // = number of columns of M
	for (col = 0; col < nCols; col++)
	    {
		if (M[col].length != nRows )
		    {   throw ("inconsistent row length at row "+ row + "of M")
			return 
		    }
	}
	result = Array(nRows)
	for (row = 0; row < nRows ; row ++)
	{
	    result[row]=Array(nCols)
	    for (col = 0; col < nCols; col++)
	    {
		result[row][col] = M[col][row]
	    }
	}
	return result
}


function matrixProduct(M1, M2)
{
	nRows = M1.length
	if (0 == nRows)
	{	return []
	}
	// else
	len = M1[0].length
	for (row = 1; row < nRows ; row ++)
	{
	    if (M1[row].length != len )
	    {   throw ("inconsistent row length at row " + row + "of M1")
		return 
	    }
	}
	if (M2.length != len)
	{	throw ("incompatible matrix sizes")
		return [] // error
	}
	// else
	nCols = M2[0].length
	for (i = 1; i<len; i ++)
	{	if (M2[i].length != nCols)	
	    {   throw ("inconsistent row length at row " +i+ "of M2")
		return 
	    }
	}

	var result = Array(nRows)
	for (row = 0; row < nRows ; row ++)
	{
	    result[row]=Array(nCols).fill(0)
	    for (col = 0; col < nCols; col ++)
	    {   
		for (let i = 0 ; i < len ; i++)
		{	result[row][col] += M1 [row][ i] * M2 [i][col]
		}
	    }
 	}
	return result
}

function diagonalMatrixProduct(diagonal, M2)
{
	nRows = diagonal.length
	if (0 == nRows)
	{	return []
	}
	// else
	if (M2.length != nRows)	
	    {   throw ("inconsistent number of rows of M2: "+M2.length+". diagonal: "+diagonal.length)
		return 
	    }
	nCols = M2[0].length
	if (0 == nCols)
	{	return new Array(nRows).fill(0).map(() => new Array(nCols))
	}
	for (i = 1; i<M2.length; i ++)
	{
	    if (M2[i].length != nCols)	
	    {   throw ("inconsistent row length at row "+ i + "of M2")
		return 
	    }
	}
	
	var result = Array(nCols)
	for (row = 0; row < nRows ; row ++)
	{
	    result[row]=Array(nCols).fill(0)
	    for (col = 0; col < nCols; col ++)
	    {   
		{	result[row][col] = diagonal [row] * M2 [row][col]
		}
	    }
 	}
	return result
}

function innerProduct(v1, v2, IPmetric)
{
	v1T = matrixProduct([v1],IPmetric)
	return matrixProduct([v2],v1T)[0][0]
}

function innerProductForDiagonalIP(v1, v2, diagonal)
{
	if ((v1.length != v2.length) || (v1.length != diagonal.length))
	{	throw ("incompatible vector lengths")
	}
	var result = 0
	for (i =0; i<v2.length; i++)
	{	result += v1[i]*v2[i]*diagonal[i]
	}
	return result
}

function innerProductForUnitIP(v1,v2)
{
	if (v1.length != v2.length)
	{	throw ("incompatible vector lengths")
	}
	var result = 0
	for (i =0; i<v2.length; i++)
	{	result += v1[i]*v2[i]
	}
	return result
}
function leastSquares(matrix, vector, IPmetric) {
  return leastSquaresForUnitIP(matrix, vector, matrixProduct(IPmetric, matrix))
}

function leastSquaresForDiagonalIP(matrix, vector, diagonalElements) {
  return leastSquaresForUnitIP(matrix, vector, 
		diagonalMatrixProduct(diagonalElements, matrix))
}

// Function to compute the least-squares solution to an overdetermined system.
// If metricA == null, ATA is calculated as transpose(A) * A
// Otherwise, ATA will be calculated as transpose(metricA) * A
function leastSquaresForUnitIP(A, B, metricA = null) {
  const m = A.length;  // Number of rows
  const n = A[0].length;  // Number of columns

  if (null == metricA)
  {	metricA = A
  }
  // Step 1: Compute A^T * A (the normal matrix)
  let ATA = Array.from(Array(n), () => Array(n).fill(0));
  // use map to make sure that not all rows will refer to the same elements.
  // let ATA = new Array(n).map(() => new Array(n).fill(0));

  for (let i = 0; i < n; i++) {
    for (let j = 0; j < n; j++) {
      for (let k = 0; k < m; k++) {
        ATA[i][j] += metricA[k][i] * A[k][j];
      }
    }
  }

  // Step 2: Compute A^T * B
  let ATB = new Array(n).fill(0);
  for (let i = 0; i < n; i++) {
    for (let k = 0; k < m; k++) {
      ATB[i] += metricA[k][i] * B[k];
    }
  }

  // Step 3: Solve the normal equation (A^T * A) * X = A^T * B
  return solveMatrixEquation(ATA, ATB);
}


// Example usage
let A = [
  [2, 1],
  [1, 3],
  [3, 2],
  [4, 1]
];

let B = [5, 6, 7, 8];

const leastSquaresSolution = leastSquaresForUnitIP(A, B);
console.log("========== Least Squares Solution:", leastSquaresSolution);

console.log("========== matrix with row and column pivoting:")
// matrix designed for row and column pivoting
A = [
  [2, 4, 4, 7, 0],
  [1, 6, 3, 2, 1],
  [2, 9, 0, 1, 4],
  [3, 2, 5, 0, 8],
  [4, 1, 3, 3, 9]
];
findBasis(A)
B = [5, 6, 7, 8, 9];

console.log("correct solution: "+solveMatrixEquation(A,B))
X = solveMatrixEquation(A,B)
console.log("experimental solution: "+X)
console.log("controle: "+B+" = "+X.multiplyByMatrix(A))

console.log("========== now a singular matrix:")
// now a singular matrix - with a possible solution
A = [
  [2, 4, 4, 7, 0],
  [3, 2, 5, 0, 8],
  [1, 6, 3, 2, 1],
  [2, 9, 0, 1, 4],
  [4, 8, 8, 2, 9],
];
findBasis(A)
B = [5, 6, 7, 8, 9].multiplyByMatrix(A)
console.log("B for square singular matrix: "+B)
X = solveMatrixEquation(A,B)
console.log(" X ="+ X)
console.log("controle: "+B+" = "+X.multiplyByMatrix(A))

// now a singular matrix - with a possible solution

B = [5, 6, 7, 8, 9].multiplyByMatrix(A)
console.log("========== B for 6 by 5 rectangular singular matrix of rank 4: "+B)
A = [
  [2, 4, 4, 7, 0],
  [3, 2, 5, 0, 8],
  [1, 6, 3, 2, 1],
  [2, 9, 0, 1, 4],
  [2, 9, 0, 1, 4],
  [4, 8, 8, 2, 9],
];
findBasis(A)
X = solveMatrixEquation(A,B)
console.log("X ="+ X)
console.log("controle: "+B+" = "+X.multiplyByMatrix(A))

console.log("==========This singular 6x5 matrix equation should have a solution:")
// now a singular of rank = nColumns that has a solution
A = [
  [2, 4, 4, 7, 0],
  [1, 6, 3, 2, 1],
  [2, 9, 0, 1, 4],
  [3, 2, 5, 0, 8],
  [4, 1, 3, 3, 9],
  [3, 1, 3, 5, 2]
];
findBasis(A)
B = [5, 6, 7, 8, 9].multiplyByMatrix(A);
X = solveMatrixEquation(A,B)
console.log("X ="+ X)
console.log("controle: "+B+" = "+X.multiplyByMatrix(A))

B = [118,87,108,134,153,100];
console.log("This  equation should not have a solution:")
X = solveMatrixEquation(A,B)
console.log("X ="+ X)
console.log("controle: "+B+" = "+X.multiplyByMatrix(A))

/*
*/
