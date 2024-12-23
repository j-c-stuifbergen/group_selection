// requires
// <script src="./js/matrix.js"></script>

function download(filename, text) {
  var element = document.createElement('a');
  element.setAttribute('href', 'data:text/plain;charset=utf-8,' + encodeURIComponent(text));
  element.setAttribute('download', filename);

  element.style.display = 'none';
  document.body.appendChild(element);

  element.click();

  document.body.removeChild(element);
}

function shuffle(array) 
{  let currentIndex = array.length;

  // While there remain elements to shuffle...
  while (currentIndex != 0) {

    // Pick a remaining element...
    let randomIndex = Math.floor(Math.random() * currentIndex);
    currentIndex--;

    // And swap it with the current element.
        [array[currentIndex], array[randomIndex]] = [
      array[randomIndex], array[currentIndex]];
  }
  
	return array
}
function sumArray(q)
{	result = 0
	for (let i=0; i<q.length ; i++)
	{	result += q[i]
	}
	return result
}

function htmlFromArray(q)
{	return "[" + q + "]" 
}

function groups(nPlaces, groupInfo) // groupSizes, nGroups) 
{
	this.nPlaces = nPlaces  // e.g. 80
	this.groupSizes = groupInfo.map(item => item.groupSize) // e.g. [1, 2, 3, 5] : groups can be 1, 2, 3 persons
	this.nGroups = groupInfo.map(item => item.nGroups)  // e.g. [40, 8, 2, 1] for 40 individuals, 8 duos, 2 trios...
	this.nCandidates = 0
	this.nIndividuals = [] // groupSizes * nGroups
	this.indivProba = [] // individual probability
	this.indivProbaDevi = [] // relative deviation from the average probability
	this.combinations = [] // an array of combinations that fill all places
	this.probabilities = [] // same dimensions as this.combinations. 
	this.differences = [] // differences with proba of first selected combination
	this.selection = [] // can contain indices that refer to this.combinations
	this.availabel = [] // the indices that are not in selection
	this.p = [] // weights, belonging to the selection.
	this.setIndividuals() // also fills this.nCandidates
	this.pAim = this.nPlaces/this.nCandidates; //
}

groups.prototype.setProbabilities = function()
// calculate the probabilities for groups in a combination
{	this.probabilities = Array(this.combinations.length)
	for (let i = 0 ; i < this.combinations.length; i++)
	{	this.probabilities[i] = Array(this.groupSizes.length) 
		for (let j = 0 ; j< this.groupSizes.length ; j++)
		{
			this.probabilities[i][j] = 
			this.combinations[i][j] / this.nGroups[j]
		}
	}
}

groups.prototype.makeDistribution = function()
{
/* make the total of this.p[] = 1
 * required: 0<=this.p[i] for each i
 * required: 0<sum of this.p[i]
 */
	var coefficientsAreValid = true
	var pTotal =0;
	for (let i = 0; i< this.p.length ; i++)
	{	if(this.p[i] <0)
		{	this.p[i]=0
			 coefficientsAreValid = false
		}
		pTotal+=this.p[i];
	}	
	for (let i = 0; i< this.p.length ; i++)
	{	this.p[i]=this.p[i]/pTotal;
	}
	return coefficientsAreValid
}

groups.prototype.setIndividuals = function()
{	this.nIndividuals = new Array(this.groupSizes.length)

	this.nCandidates = 0
	for (let i=0; i<this.groupSizes.length ; i++)
	{	this.nIndividuals[i] = this.groupSizes[i] * this.nGroups[i]
		this.nCandidates += this.nIndividuals[i]
	}
	return this.nIndividuals
}

groups.prototype.parabolaMin = function(fie, x=1.0)
{
	// y = c0 + c1 x + c2 x^2
	delta = 0.5*x
	y0 = fie(x-delta)
	y1 = fie(x)
	y2 = fie(x+delta)

	c2 = (y0-2*y1+y2) / (4*delta)

	y1 = y1-c2*x*x
	c1 = (y1 - (y0-c2*(x-delta)*(x-delta))) / delta

	//c0 = y1-c1*x

	return c1 / (-2.0*c2)
}


/*
	// definitions
	// nIndividuals[i] = nGroups[i]*groupSizes[i] 
	// q[n][i] = v[n] [i] / nGroups[i] 
	// p[i] = sum_k  w[k] q[k][i]  
	// f = sum_j    ( p[j] - a[j] )^2  * nIndividuals[j]
	
	 first derivatives
	// dp[i] / dw[k] = q[k][i]
	// df / dp[i] = 2 (p[i]-a[i]) * nIndividuals[i]
			= 2 ( sum_k  w[k] q[k][i]  - a[i] ) *nIndividuals[i]
	// df / dw[n] = sum_i ( df / dp[i] * dp[i]/dw[n] )
		= sum_i ( 2  (sum_k  w[k] q[k][i] - a[i] ) *nIndividuals[i]
		  * q[n][i] )
     	// if 0 == df / dw[n] then for every n:
      	sum_i (sum_k w[k] q[k][i] * q[n][i] * nIndividuals[i] = sum_i a[i] * q[n][i] * nIndividuals[i]
       This is a matrix equation for w:
       Q w = b
       Q[k,n] = < q[k], q[n] > (inner product defined by nIndividuals[])
       b[n]   = < a   , q[n] >
       
	// d df / dw[n] dw[m] = sum_i ( 2   q[m][i]   *nIndividuals[i]
		  * q[n][i] )
	// 
*/	
groups.prototype.findMinimum = function()
/* find coefficients that minimize the difference between realized and aimed probabilities,
without requiring that the probabilities will add up to 1
Note: such a minimum always exists.
/* This is a matrix equation for the weights:
       Q p= b
       Q[k,n] = < q[k], q[n] > (inner product defined by nIndividuals[])
       b[n]   = < a   , q[n] >
 Note that the matrix can be singular.   
*/
{	var aimVector = new Array(this.groupSizes.length).fill(this.pAim)
	var b = Array(this.selection.length)
	var Q=Array(this.selection.length).fill(0).map(() => new Array(this.selection.length))
	for (var n=0; n< this.selection.length ; n++)
	{	b[n] = this.innerProd (aimVector, this.probabilities[this.selection[n]])
	 	for(var k=n; k<this.selection.length; k++)
			{	Q[k][n] = this.innerProd (this.probabilities[this.selection[k]], 
							     	this.probabilities[this.selection[n]])
				Q[n][k] = Q[k][n]
			}
	}
	this.p = solveMatrixEquation(Q,b)
	/*console.log("Q is \n")
	for (let i =0; i<Q.length; i++)
	console.log(Q[i])
	console.log("b is \n"+b)

	console.log("solution is \n"+this.p)
	*/
}
groups.prototype.testExistence = function()
/* find coefficients that minimize the difference between realized and aimed probabilities,
without requiring that the probabilities will add up to 1
Note: such a minimum always exists.
We require that the probabilities add up to 1, so we start with the best vector ("b")
We then subtract this from the aimed probability ("a"), and perform a Least Squares routine
to calculate the nearest solution.

   This is a matrix equation for the weights:
       Q (p-b)= c
       Q[k,n] = < q[k]-b, q[n]-b > (inner product defined by nIndividuals[])
       c[n]   = < a-b   , q[n]-b >
 Note that the matrix can be singular, but still has a solution.   
*/
{	var aimVector = new Array(this.groupSizes.length).fill(this.pAim)

	this.availabel = Array.from(Array(this.combinations.length).keys())

	this.selectLowestPenalty(true)
	this.makeDifferences()
	var aimDifference = aimVector.subtract(this.probabilities[this.selection[0]])

	var c = Array(this.differences.length)
	var Q=Array(this.differences.length).fill(0).map(() => new Array(this.differences.length))
	for (var n=0; n< this.differences.length ; n++)
	{	c[n] = this.innerProd (aimDifference, this.differences[n])
	 	for(var k=n; k<this.differences.length; k++)
			{	Q[k][n] = this.innerProd (this.differences[k], 
							     	this.differences[n])
				Q[n][k] = Q[k][n]
			}
	}
	pPlus = solveMatrixEquation(Q,c)
	console.log("pPlus is "+pPlus)
	console.log("pPlus.length = "+pPlus.length)
	console.log("probaOfBest + aimDifference = "+this.probabilities[this.selection[0]]+" + "+aimDifference+" = "+aimVector+" = aimVector")	
	pBest = 1
	for (let j = 0; j < pPlus.length ; j++)
	{	pBest -= pPlus[j]
	}
	
	this.p = Array(this.combinations.length).fill(0)
	for (let j=0; j < this.availabel.length ; j++)
	{	this.p[this.availabel[j]]= pPlus[j]
		this.selection.push(this.availabel[j])
	}
	this.p[this.selection[0]] = pBest

	pPlus.unshift(pBest)
	console.log("pPlus na unshfit is "+pPlus)
	console.log("pPlus.length = "+pPlus.length)

	console.log("solution if negative probabilities were allowed: "+pPlus)
	this.p = pPlus
}

groups.prototype.makeDifferences = function(probaVector= null, selection = null)
{
	if (null == probaVector)
	{	probaVector = this.probabilities[this.selection[0]]
	}
	if (null == selection)
	{	selection = this.availabel
	}
	// console.log("differences")
	this.differences = Array(selection.length)
	for (let j = 0 ; j< this.differences.length; j++)
	{	this.differences[j] = this.probabilities[selection[j]].subtract(probaVector)
		// console.log(this.differences[j])
	}
}

groups.prototype.penalty = function(combi)
{	var result = 0
	for(let i = 0 ; i < this.groupSizes.length ; i++)
	{	// delta^2 * number of individuals.
		result += Math.pow(this.pAim - combi[i], 2) * this.nGroups[i] * this.groupSizes[i] 
	}
	return result
}

groups.prototype.derivative = function()
/*	// df / dw[n] = sum_i ( df / dp[i] * dp[i]/dw[n] )
		= sum_i ( 2  (sum_k  w[k] q[k][i] - a[i] ) *nIndividuals[i]
		  * q[n][i] )
*/
{	var q = []; // used to adapt p, store deviations from desired probabilities
	// calculate probabilities per individual 
	let p = this.pIndividuals() // p[i] = sum_k (  w[k] q[k][i])
	console.log("kansen: "+ p)

	// calculate difference from desired probability
	for (let i = 0; i< this.selection.length ; i++)
	{	q[i]=0;
		for (let j = 0; j<this.groupSizes.length; j++)
		{
		  // give extra weight if a deviation affects more persons
		  q[i]+=2 * (this.pAim - p[j])*this.combinations[this.selection[i]][j]*this.groupSizes[j];
		  // don't make chances negative
		  if (q[i] < 0 && 0 ===this.p[i])
		  {	q[i] = 0
		  }
		}
	}
	return q
}

groups.prototype.innerProd = function(prob1, prob2)
// inner product , weighted by Groupsizes*nGroups
{	
	return innerProductForDiagonalIP(prob1, prob2, this.nIndividuals)
}

groups.prototype.selectNSimilarCombinations = function(N = 6)
// this function tends to select combinations that are too similar.
{	this.selection = [0]
	
	penaltyArr = [this.penalty(this.probabilities[0])]
	indexOfWorst = 0
	
	for (var i =1; i<N ; i++)
	{
		this.selection[i] = i
		penaltyArr[i] = this.penalty(this.probabilities[i])
		if (penaltyArr[indexOfWorst]<penaltyArr[i]  )
		{	indexOfWorst = i
		}
	}
	for(; i<this.combinations.length ; i++)
	{	penalty = this.penalty(this.probabilities[i])
		if (penalty < penaltyArr[indexOfWorst])
		{	
			this.selection[indexOfWorst] = i
			penaltyArr[indexOfWorst] = penalty

			// find the next worst...
			for (j = 0; j<penaltyArr.length; j++)
			{	if (penalty < penaltyArr[j])
				{ 	indexOfWorst = j
					penalty = penaltyArr[j]
				}
			}	
		}
	}
}

groups.prototype.selectBestCombinations = function(epsilon = 1e-6)
{	console.log("=== find the best combinations by least squares ===")

	aimVector = new Array(this.groupSizes.length).fill(this.pAim)
	console.log("aimVector " +aimVector)
	this.setProbabilities()

	pIndividu = Array(this.groupSizes.length).fill(0)

	this.selection = []
	this.availabel = Array.from(Array(this.combinations.length).keys())
	this.p = []		
	
	this.selectLowestPenalty(true) 
	this.makeDifferences()
	var aimDifference = aimVector.subtract(this.probabilities[this.selection[0]])
	// perform projection, so the inner product with pDiff will be 0
	this.p = leastSquaresForDiagonalIP(
		this.probabilityMatrix(), aimVector, this.nIndividuals )

	while (true)
	{	 
		pIndividu.fill(0)
		for (let i=0; i<this.p.length ; i++)	
		{	pIndividu = pIndividu.add(
			this.probabilities[ this.selection[i]]. timesScalar(this.p[i]))
		}

		pDiff = aimVector.subtract(pIndividu)
		penalty = this.innerProd(pDiff, pDiff)
		// penalty = this.penalty(pIndividu)
		console.log("weights for the combinations are " + this.p)
		console.log("penalty =" + penalty + " ----- pDiff =" + pDiff)
		
		if ( penalty < epsilon) // objectives reached!
		{  	break
		}
		if (this.p.length >= this.groupSizes.length) 
		{ 	// I can't add another vector for Least Squares
			break
		}

		ipResult = this.selectByIP(pDiff, true)
		if (ipResult.index < 0)
		{	// no vector has been found
			break
		}
		if (epsilon < ipResult.product )
		{	this.p = leastSquaresForDiagonalIP(
				this.probabilityMatrix(), aimVector, this.nIndividuals )
		}
		else
		{	// The new vector doesn't contribute significantly
			this.p.push(0)
			break
		}
	}
	let noNegatives = this.makeDistribution() 		
	return noNegatives && (penalty < epsilon)
}

groups.prototype.selectBestCombinationsOld = function(epsilon = 1e-6)
{	console.log("=== find the best combinations by least squares ===")

	aimVector = new Array(this.groupSizes.length).fill(this.pAim)
	console.log("aimVector " +aimVector)
	this.setProbabilities()

	pIndividu = Array(this.groupSizes.length).fill(0)

	this.selection = []
	this.availabel = Array.from(Array(this.combinations.length).keys())
	this.p = []		
	
	if (true)  // perhaps not necessary?
	{ this.selectLowestPenalty(true) 
	// perform projection, so the inner product with pDiff will be 0
	this.p = leastSquaresForDiagonalIP(
		this.probabilityMatrix(), aimVector, this.nIndividuals )
	}
	while (true)
	{	 
		pIndividu.fill(0)
		for (let i=0; i<this.p.length ; i++)	
		{	pIndividu = pIndividu.add(
			this.probabilities[ this.selection[i]]. timesScalar(this.p[i]))
		}

		pDiff = aimVector.subtract(pIndividu)
		penalty = this.innerProd(pDiff, pDiff)
		// penalty = this.penalty(pIndividu)
		console.log("weights for the combinations are " + this.p)
		console.log("penalty =" + penalty + " ----- pDiff =" + pDiff)
		
		if ( penalty < epsilon) // objectives reached!
		{  	break
		}
		if (this.p.length >= this.groupSizes.length) 
		{ 	// I can't add another vector for Least Squares
			break
		}

		ipResult = this.selectByIP(pDiff, true)
		if (ipResult.index < 0)
		{	// no vector has been found
			break
		}
		if (epsilon < ipResult.product )
		{	this.p = leastSquaresForDiagonalIP(
				this.probabilityMatrix(), aimVector, this.nIndividuals )
		}
		else
		{	// The new vector doesn't contribute significantly
			this.p.push(0)
			break
		}
	}
	let noNegatives = this.makeDistribution() 		
	return noNegatives && (penalty < epsilon)
}

groups.prototype.selectLowestPenalty = function(update = false )
{	
	if (0 == this.availabel.length)
	{	return { penalty: null, index:-1}
	}

	penalty = this.penalty(this.probabilities[this.availabel[0]])
	indexOfLowestPenalty = this.availabel[0];  
	var foundAt

	for(let i = 1; i<this.availabel.length ; i++)
	{	newPenalty = this.penalty(this.probabilities[this.availabel[i]])
		if (newPenalty < penalty)
		{	penalty = newPenalty
			indexOfLowestPenalty =this.availabel[i] 
			foundAt = i
		}
	}
	if (update && ( 0 <= indexOfLowestPenalty))
	{	this.selection = [indexOfLowestPenalty]
		this.availabel.splice(foundAt,1)
	}
	return {penalty: penalty, index: indexOfLowestPenalty}
}

groups.prototype.selectByIP = function(goalVector, update = false)
{		 
	indexOfBestCombination = -1; // 
	biggestIP = 0 ;

	for (let i=0; i<this.availabel.length ; i++)
	{	ip = this.innerProd(goalVector, this.probabilities[this.availabel[i]])
		if (biggestIP < ip)
		{	biggestIP = ip
			indexOfBestCombination = this.availabel[i]
			foundAt     = i
		}
	}
	if (update && ( 0 <= indexOfBestCombination))
	{	this.selection.push(indexOfBestCombination)
		this.availabel.splice(foundAt, 1)
	}
	return {product: biggestIP, index:indexOfBestCombination}
}

groups.prototype.solveLeastSquares = function(aimVector)
{
	this.p = leastSquaresForDiagonalIP(
		this.probabilityMatrix(), aimVector, this.nCandidates )

	pIndividu.fill(0)
	for (let i=0; i<this.p.length ; i++)	
	{	pIndividu.add(
		this.probabilities[this.selection[i]].timesScalar(this.p[i]))
	}
	return this.p
}

groups.prototype.initialWeights = function()
{	
	this.p=[];
	for (let i = 0; i< this.selection.length ; i++)
	{	this.p[i]=1/this.selection.length;
	}
	
	return this.p
}
groups.prototype.improveWeights = function(nIterations, stepFactor = 1.001, htmlId = null)
// at each iteration, adaptfactor will be multiplied by stepfactor
{
	let validity = true
	if (null == this.p)
	{	this.p = this.initialWeights()
	}
	var maxWeight = 0;
	for(let i = 0; i<this.groupSizes.length; i++)
	{	maxWeight = Math.max(maxWeight,this.groupSizes[i]*this.nGroups[i])
	}

	var adaptFactor = 0.5/maxWeight/this.selection.length
	console.log("adaptFactor = "+adaptFactor)
		
	for (let k = 0; k< nIterations ; k++)
	{
		q= this.derivative()
		q2 = 0
		for(let i=0; i<q.length; i++)
		{	q2 += q[i]*q[i];
		}
		message = "iteration "+k+", the norm of the derivative is "+Math.sqrt(q2) +"<br>\n"
		var pTotal = 0
		
		validity = true
		for (let i = 0; i< this.selection.length ; i++)
		{	this.p[i]+=(q[i])*adaptFactor
			
			// correct for overshoot
			if (this.p[i]<0)
			{	this.p[i]=0
				validity = false
			}
			pTotal+=this.p[i]
		}
		for (let i = 0; i< this.selection.length ; i++)
		{	
			this.p[i]=this.p[i]/pTotal
		}
		adaptFactor *= stepFactor
	}
	return validity
}

groups.prototype.improveWeightsOud = function(nIterations, stepFactor = 1.001)
// at each iteration, adaptfactor will be multiplied by stepfactor
{
	if (null == this.p)
	{	this.p = this.initialWeights()
	}
	var maxWeight = 0;
	for(let i = 0; i<this.groupSizes.length; i++)
	{	maxWeight = Math.max(maxWeight,this.groupSizes[i]*this.nGroups[i])
	}

	var adaptFactor = 0.5/maxWeight/this.selection.length
	console.log("adaptFactor = "+adaptFactor)
		
	for (let k = 0; k< nIterations ; k++)
	{
		q= this.derivative()
		q2 = 0
		for(let i=0; i<q.length; i++)
		{	q2 += q[i]*q[i];
		}
		console.log("iteration "+k+", the norm of the derivative is "+Math.sqrt(q2))
		var pTotal = 0
		for (let i = 0; i< this.selection.length ; i++)
		{	this.p[i]+=(q[i])*adaptFactor
			
			// correct for overshoot
			if (this.p[i]<0)
			{	this.p[i]=0
			}
			pTotal+=this.p[i]
		}
		console.log("coeff : "+this.p)
		for (let i = 0; i< this.selection.length ; i++)
		{	
			// this.p[i]=this.p[i]/pTotal
		}
		// p = this.makePdistribution(p)
		adaptFactor *= stepFactor
	}
	return this.p
}

groups.prototype.selectionArray= function(selection = this.selection, combis=this.combinations)
{
	result = Array(selection.length)
	for (let i = 0; i< selection.length; i++)
	{	result [i] = combis[selection[i]]
	}
	return result
}

groups.prototype.probabilityMatrix= function()
{	return  transpose(this.probabilityArray())
}

groups.prototype.probabilityArray= function()
{	return this.selectionArray(this.selection, this.probabilities)
}

groups.prototype.pIndividuals = function()
{	// p contains the probability of vectors in de selection
	var e = [];// expected number of groups that will have a place
	this.indivProba = []; // p[j]=e[j]/nGroups[j]; probability that a member of these groups will have a place
	for (let j = 0; j<this.groupSizes.length; j++)
	{
		e[j] = 0; // expected number of groups of groupSize[j]
		for (let i = 0; i< this.selection.length ; i++)
		{	e[j] += this.p[i]*this.combinations[this.selection[i]][j];
		}
		this.indivProba[j]=e[j]/this.nGroups[j];
		this.indivProbaDevi[j]=this.indivProba[j]/this.pAim -1
	}
	return this.indivProba
}
groups.prototype.selectAllCombinations = function()
{	this.selection = Array.from(Array(this.combinations.length).keys())
}

groups.prototype.setCombinations = function()
{
	const recursivePart = (nPlaces, nGroupsIn) =>
	{
		var nGroups = []
		for(let i=0;i<nGroupsIn.length;i++)
		{	nGroups.push(nGroupsIn[i])
		}
		var max = nGroups.pop()
		if(0==nGroups.length)
		{
			if ( max*this.groupSizes[0] < nPlaces)
			{	// cannot fill all places
				return []
			}
			if (0==nPlaces%this.groupSizes[0])
			{
				return [[nPlaces/this.groupSizes[0]]]
			}
			// cannot fill all places without splitting a group
			return []
		}
		// else
		var result = []
			
		for (let i = 0; i<=max; i++)
		{	nLeft = nPlaces - i*this.groupSizes[nGroups.length]
			if (0<=nLeft)
			{
				subset = recursivePart(nLeft, nGroups)
				for (let j=0;j<subset.length;j++)
				{	subset[j].push(i)
				}
				result = result.concat(subset)
			}
		}
		return result
	}

	this.combinations = recursivePart(this.nPlaces,this.nGroups)
}
groups.prototype.calculationString = function(index, selec = this.selection)
{	
	if (selec.length <= index)
	{	return "There are not more than " + selec.length +" selected combinations."
			+ i + " is too big."
	}
	combi = this.combinations[selec[index]]
	var result = ""
	let total = 0
	let j = 0
	for (; j<this.groupSizes.length-1; j++)
	{	total +=  combi[j] * this.groupSizes[j]
		result += combi[j]+"x"+this.groupSizes[j]+" + "
	}
	total +=  combi[j] * this.groupSizes[j]
	result += combi[j]+"x"+this.groupSizes[j]+" = " + total
	
	return result 
}

groups.prototype.calculationGroupsString = function()
{	
	return this.calculationVectorToTeX ( this.nGroups)
}
groups.prototype.selectedToTex = function (index)
{
	return vector_to_TeX (this.combinations[this.selection[index]])
}
groups.prototype.calculationTeX = function(i, includeTot = true, selec = this.selection)
{	
	if (selec.length <= i)
	{	return "There are not more than " + selec.length +" selected combinations."
			+ i + " is too big."
	}
	return this.calculationVectorToTeX ( this.combinations[selec[i]], includeTot)
}

groups.prototype.calculationVectorToTeX = function(combi, includeTotal = true)
{	
	var result = "\\begin{pmatrix}"
	let total = 0
	for (j=0; j<this.groupSizes.length; j++)
	{	total +=  combi[j] * this.groupSizes[j]
		result += 0<j? "\\\\ " : ""
		result +=  "\\text{"+combi[j]+"x"+this.groupSizes[j] + "} &=& "
				+ combi[j] * this.groupSizes[j]
	}
	if(includeTotal)
	{	result += "\\\\ \\text{total} &=& " + total 
	}
	result += "\\end{pmatrix} "
	
	return result 
}

groups.prototype.chancesTable = function(nDigits = 3)
{	
	this.pIndividuals()
	
	var result = "<table>"
	result +=    "<tr><th>Groups of size</th>"
		+    "<th>probability for persons<br>in this group </th>"
		+    "<th>rel. error = <br>proba/average - 1</th>"
		+    "<th>expected number <br> of groups of this size</th>"
	result +=    "<th>Expected number <br>of persons</th></tr>"

	var totalCandidates = 0
	var totalAdmitted = 0
	var rtsqError = 0
	for (j=0; j<this.groupSizes.length; j++)
	{	let nCandidates = (this.nGroups[j]*this.groupSizes[j])
		let expectedGroups = (this.indivProba[j]*this.nGroups[j])
		let expectedAdmitted = expectedGroups*this.groupSizes[j]
		totalCandidates += nCandidates
		totalAdmitted += (expectedAdmitted)
		result += "<tr><td>" + this.nGroups[j]+" x " + this.groupSizes[j] 
			+ " person(s) = " + nCandidates +"</td>"
			+ "<td align='right'>"+this.indivProba[j].toFixed(8) +"</td>"
			+ "<td align='right'>"+this.indivProbaDevi[j].toFixed(nDigits) +"</td>"
			+ "<td align='right'>" 
			+ this.nGroups[j] +	" x "+this.indivProba[j].toFixed(nDigits) 
			+ " = " 	+ expectedGroups.toFixed(nDigits)
			+ "</td><td align='right'>" +expectedGroups.toFixed(nDigits) + " x " 
			+ (this.groupSizes[j]) + " = "
			+ expectedAdmitted.toFixed(nDigits)
			+ "</td> </tr> \n"
		rtsqError += Math.pow(this.indivProba[j]-this.pAim,2) * this.nIndividuals[j]
	}
	result += "<tr><th align='right'> total: "+this.nCandidates+"</th>" + //+totalCandidates+"</th>" +
		  "<th align='right'>average P: "+this.pAim.toFixed(nDigits) + "</th>" +
		  "<th align='right'>weighted error: "+rtsqError.toFixed(nDigits) + "</th>" +
		"<th></th>" +
		  "<th align='right'> total: "+totalAdmitted.toFixed(9) + "</th>" +
		"</table>"
	
	return result 
}

groups.prototype.resultsTable = function(nDigits = 8)
{	
	var result = "<table>"
	result +=    "<tr>"
                + 	"<th>number of<br>persons</th>"
                + 	"<th>Groups of size</th>"
                + 	"<th>number of groups</th>"
	for (i=0; i<this.selection.length; i++)
	{	result += "<th>combination "+i+"<br>probability:<br>"
		+this.p[i].toFixed(nDigits)+"</th>"
	}
	result +="</tr>\n"
	for (j =0; j<this.groupSizes.length; j++)
	{
		result +=    "<tr><td align='right'>" + (this.groupSizes[j]*this.nGroups[j]) +"</td> "
		result +=    "<td align='right'>" + (this.groupSizes[j]) +"</td> "
		result +=    "<td align='right'>" + (this.nGroups[j]) +"</td> "
		for (i=0; i<this.selection.length; i++)
		{	result += "<td align='right'>"
			 + this.combinations[this.selection[i]][j] +"</td>"
		}
		result += "</tr>\n"
	}
	result +=    "<tr>"
                + 	"<td> </td>"
                + 	"<td> </td>"
                + 	"<th>probability</th>"
	for (i=0; i<this.selection.length; i++)
	{	result += "<th>"+this.p[i]+"</th>"
	}
	result +="</tr>\n"
	result +=    "<tr>"
                + 	"<td> </td>"
                + 	"<td> </td>"
                + 	"<th>cumulative(sum)</th>"
	let total = 0
	for (i=0; i<this.selection.length; i++)
	{	total += this.p[i]
		result += "<th>"+total+"</th>"
	}
	result +="</tr>\n"

	result += 	"</table>\n"
	
	return result 
}

groups.prototype.resultCsv = function(selec = this.selection, p=this.p, separator = "\t", newline = "\r\n")
{
	result = ""
	// add headers
	result += this.nPlaces + separator + "'= total places "
		+ separator + "combination index-->"
	for (let i = 0; i<p.length; i++)
	{	result += separator + i
	}	
	// add Weights
	result += newline + this.nCandidates + separator + "'=total persons " 
			+ separator + "probability of combi-->"
	for (let i = 0; i<p.length; i++)
	{	result += separator + p[i]
	}
	// add separator
	result += newline+ "------ " + separator + "------ " 
		+ separator + "------ "
	for (let i = 0; i<p.length; i++)
	{	result += separator + "-----"
	}
	// add group header   
	result += newline+"N_persons " + separator + "groupsize " +
		 separator + "N_groups" + separator + "selection -->"
	// add combinations
	for (j = 0; j<this.groupSizes.length; j++)
	{   result += newline + this.nIndividuals[j]
		 + separator + this.groupSizes[j] + separator + this.nGroups[j] 
	    for (let i = 0; i<p.length; i++)
	    {	result += separator + this.combinations[selec[i]][j]
	    }
	}
	// optional: add sums?
	result += newline
	result += newline + "individual probability" + separator +
			"expected persons"+ separator +"expected groups " 
			
	colNames = Array.from("ABCDEFGHIJKLMNOPQRSTU")
	endColumn = colNames[2+this.selection.length]
	startColIndex = 3
	startRow = 7+this.groupSizes.length
	for (j = 0; j<this.groupSizes.length; j++)
	{   result += newline + "=$C"+(startRow+j)+"/$C"+(5+j)
		+ separator +  "=$C"+(startRow+j)+"*$B"+(5+j) + separator 
			+ "=SUM($D"+(startRow+j)+":$"+endColumn+(startRow+j)+")"
			
	    for (let i = 0; i<p.length; i++)
	    {	result += separator + "="+colNames[startColIndex + i]+"$2*" 
				+colNames[startColIndex + i]+(5+j) 
			//+this.combinations[selec[i]][j]*this.p[i]
	    }
	}
		 
	return result.replaceAll('.', ',')
}

groups.prototype.selectionString = function(selec = this.selection, newLine = "<br>") // for console: "\r\n"
{
	let combinationString = ""
	let calcStr = ""
	for (let i = 0; i<selec.length; i++)
	{	combinationString += "["+this.combinations[selec[i]]+"]"+newLine
		calcStr+= this.calculationString(i, selec)+newLine
 	}
	return {"combinations":combinationString, "calculations": calcStr}
}

