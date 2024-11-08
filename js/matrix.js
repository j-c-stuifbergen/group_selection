Array.prototype.timesScalar = function(scalar)
{	result = new Array(this.length)
	for(let i = 0; i<this.length ; i++)
	{	result[i] = this[i]*scalar
	}
	return result
}
Array.prototype.add = function(vector)
{	
	result = new Array(this.length)
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
	result = new Array(this.length)
	if (this.length != vector.length)
	{	throw ("addition of vectors of different length")
	}
	for (i = 0; i< this.length; i++)
	{	result[i] = this[i] - vector[i]
	}
	return result
}
// Function to solve a matrix equation A * X = B
function solveMatrixEquation(A, B) {
  const n = A.length;
  
  // Augment matrix A with column vector B
  let augmentedMatrix = A.map((row, i) => [...row, B[i]]);
  
  // Forward elimination with partial pivoting
  for (let col = 0; col < n; col++) {
    // Find the pivot row
    let pivotRow = col;
    for (let row = col + 1; row < n; row++) {
      if (Math.abs(augmentedMatrix[row][col]) > Math.abs(augmentedMatrix[pivotRow][col])) {
        pivotRow = row;
      }
    }

    // Swap the current row with the pivot row
    if (pivotRow !== col) {
      let temp = augmentedMatrix[col];
      augmentedMatrix[col] = augmentedMatrix[pivotRow];
      augmentedMatrix[pivotRow] = temp;
    }

    // Perform elimination to make all values below the pivot zero
    for (let row = col + 1; row < n; row++) {
      const factor = augmentedMatrix[row][col] / augmentedMatrix[col][col];
      for (let j = col; j <= n; j++) {
        augmentedMatrix[row][j] -= factor * augmentedMatrix[col][j];
      }
    }
  }

  // Back substitution
  let X = new Array(n).fill(0);
  for (let row = n - 1; row >= 0; row--) {
    let sum = augmentedMatrix[row][n]; // Right-hand side of the equation
    for (let col = row + 1; col < n; col++) {
      sum -= augmentedMatrix[row][col] * X[col];
    }
    X[row] = sum / augmentedMatrix[row][row];
  }

  return X;
}

function transpose(M)
{	// number of columns of the transpose
		console.log("M is "+M)
	nCols = M.length // = number of rows of M
	if (0 == nCols)
	{	return []
	}
	// else
	nRows = M[0].length // = number of columns of M
	for (col = 0; col < nCols; col++)
	    {
		    console.log( "nRows is "+nRows+ ", nCols is "+nCols+", col is "+ col)
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
		    if (M[col].length != nRows )
		    {   throw ("inconsistent row length at row "+ row + "of M")
			return 
		    }
		result[row][col] = M[col][row]
	    }
	}
	return result
}

function matrixProduct(M1, M2)
{
	result = []
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

	result = Array(nRows)
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
	result = []
	nRows = diagonal.length
	console.log("nRows is "+nRows)
	if (0 == nRows)
	{	return []
	}
	// else
	if (M2.length != nRows)	
	    {   throw ("inconsistent row length at row "+ i + "of M2")
		return 
	    }
	nCols = M2[0].length
	if (0 == nCols)
	{	return new Array(nRows).fill(0).map(() => new Array(nCols))
	}
	console.log("ncols is "+nCols)
	for (i = 1; i<M2.length; i ++)
	{
	    if (M2[i].length != nCols)	
	    {   throw ("inconsistent row length at row "+ i + "of M2")
		return 
	    }
	}
	
	result = Array(nCols)
	for (row = 0; row < nRows ; row ++)
	{
	    result[row]=Array(nCols).fill(0)
	    for (col = 0; col < nCols; col ++)
	    {   
		{	result[row][col] = diagonal [row] * M2 [row][col]
		}
	    }
 	}
	console.log("diagonal : "+diagonal)
	console.log("   x matrix: "+M2)
	console.log("diagonal x matrix: "+result)
	return result
}

function innerProduct(v1, v2, IPmetric)
{
	v1T = matrixproduct([v1],IPmetric)
	return matrixProduct([v2],v1T)[0][0]
}

function innerProductForDiagonalIP(v1, v2, diagonal)
{
	if ((v1.length != v2.length) || (v1.length != diagonal.length))
	{	throw ("incompatible vector lengths")
	}
	result = 0
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
	result = 0
	for (i =0; i<v2.length; i++)
	{	result += v1t[i]*v2[i]
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

  console.log("a " + A)
  console.log("B " + B)
  console.log("at " + metricA)
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
  console.log("ata " + ATA)
  console.log("atB " + ATB)
  // Step 3: Solve the normal equation (A^T * A) * X = A^T * B
  return solveMatrixEquation(ATA, ATB);
}

// Example usage
const A = [
  [2, 1],
  [1, 3],
  [3, 2],
  [4, 1]
];

const B = [5, 6, 7, 8];

const leastSquaresSolution = leastSquaresForUnitIP(A, B);
console.log("Least Squares Solution:", leastSquaresSolution);

