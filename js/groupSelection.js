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
//	console.log("randomIndex = "+ randomIndex)
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

function groups(nPlaces, groupSizes, nGroups) 
{
	this.nPlaces = nPlaces  // e.g. 80
	this.groupSizes = groupSizes // e.g. [1, 2, 3, 5] : groups can be 1, 2, 3 persons
	this.nGroups = nGroups  // e.g. [40, 8, 2, 1] for 40 individuals, 8 duos, 2 trios...
	this.nIndividuals = [] // groupSizes * nGroups
	this.nCandidates = 0
	this.combinations = [] // an array of combinations that fill all places
	this.probabilities = [] // same dimensions as this.combinations. 
	this.selection = [] // can contain indices that refer to this.combinations
	this.p = [] // weights, belonging to the selection.
	this.setIndividuals() // also fills this.nCandidates
	this.pAim = this.nPlaces/this.nCandidates; // also fills this.nCandidates
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
	var pTotal =0;
	for (let i = 0; i< this.p.length ; i++)
	{	pTotal+=p[i];
	}	
	for (let i = 0; i< this.p.length ; i++)
	{	this.p[i]=p[i]/pTotal;
	}
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

groups.prototype.diff2 = function(v1, v2)
{	result = 0
	for( let i = 0; i<this.groupSizes.length; i++)
	{	result += Math.pow(v1[i]-v2[i], 2) * groupSizes[i]
	}
	return result
}

groups.prototype.derivative = function()
{	var q = []; // used to adapt p, store deviations from desired probabilities
	// calculate probabilities per individual 
	let a = this.pIndividuals(this.p)
	console.log("kansen: "+ a)

	// calculate difference from desired probability
	for (let i = 0; i< this.selection.length ; i++)
	{	q[i]=0;
		for (let j = 0; j<this.groupSizes.length; j++)
		{
		  // give extra weight if a deviation affects more persons
		  q[i]+=(this.pAim - a[j])*this.combinations[this.selection[i]][j]*this.groupSizes[j];
		  // don't make chances negative
		  if (q[i] < 0 && 0 ===this.p[i])
		  {	q[i] = 0
		  }
		}
	}
	return q
}
groups.prototype.penalty = function(combi)
{	var result = 0
	for(let i = 0 ; i < this.groupSizes.length ; i++)
	{	// delta^2 * number of individuals.
		result += Math.pow(this.pAim - combi[i], 2) * this.nGroups[i] * this.groupSizes[i] 
	}
	return result
}
groups.prototype.innerProd = function(prob1, prob2)
// inner product , weighted by Groupsizes*nGroups
{	
/*	console.log("nindividuals is "+ this.nIndividuals)
	console.log("prob1 is "+ prob1)
	console.log("prob2 is "+ prob2)
	console.log("nindividuals is "+ this.nIndividuals.length)
	console.log("prob1 is "+ prob1.length)
	console.log("prob2 is "+ prob2.length)
*/
	return innerProductForDiagonalIP(prob1, prob2, this.nIndividuals)
}

groups.prototype.selectFirstCombination = function()
{	this.selection = [0]
	penalty = this.penalty(this.probabilities[0])

	for(let i = 1; i<this.combinations.length ; i++)
	{	newPenalty = this.penalty(this.probabilities[i])
		if (newPenalty < penalty)
		{	penalty = newPenalty
			this.selection[0] = i
		}
	}
}

groups.prototype.selectBestCombinationsOld = function()
{	this.selection = [0]
	penalty = this.penalty(this.probabilities[0])

	for(let i = 1; i<this.combinations.length ; i++)
	{	newPenalty = this.penalty(this.probabilities[i])
		if (newPenalty < penalty)
		{	penalty = newPenalty
			// this.selection[0] = i
			// add the index at the beginning, leave old results
			this.selection.splice(0,0,i) 
		}
	}
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
{	// const epsilon = 1e-6
	this.selectFirstCombination()
	this.p=[1];	
	console.log(this.selection + " is selectie")
	pVector = this.probabilities[this.selection[0]]
	pVector = this.probabilities[0]
	aimVector = new Array(this.groupSizes.length).fill(this.pAim)

	this.p = leastSquaresForDiagonalIP(
			this.probabilityMatrix(), aimVector, this.nIndividuals )
	pVector = this.pIndividuals(this.p)
		console.log("aimVector " +aimVector)
		console.log("nIndividuals " +this.nIndividuals)
		console.log("debug eerste p " +this.p)

	while (this.p.length < this.groupSizes.length)
	{	 

		pVerschil = aimVector.subtract(pVector)
		console.log("pVerschil =" + pVerschil)
		penalty = this.innerProd(pVerschil, pVerschil)
		// penalty = this.penalty(pVector)
		console.log("penalty =" + penalty)
		if ( penalty < epsilon)
		{ 	break
		}

		this.selection.push( this.selectByIP(pVerschil) )

		this.p = leastSquaresForDiagonalIP(
			this.probabilityMatrix(), aimVector, this.nIndividuals )
	console.log('--------------------------------------------')

		console.log("debug of p " +this.p)

		pVector.fill(0)
		for (let i=0; i<this.p.length ; i++)	
		{	pVector = pVector.add(
			this.probabilities[this.selection[i]].timesScalar(this.p[i]))
	console.log("i = "+i +" proba is "+		this.probabilities[this.selection[i]].timesScalar(this.p[i]))
	console.log("i = "+i +" pVector is "+		pVector)
	console.log("coefficients are "+ this.p)
		}
	console.log('====================================================')
	}
	return pVector
}

groups.prototype.selectByIP = function(goalVector)
{		 
	bestIndex = 0; bestIP = 0
	let i = 0;
	for (; i<this.combinations.length ; i++)
	{	ip = this.innerProd(goalVector, this.probabilities[i])
		if (bestIP < ip)
		{	bestIP = ip
			bestIndex = i
		}
	}
	console.log("selected index by IP: "+bestIndex)
	return bestIndex
}

groups.prototype.solveLeastSquares = function(aimVector)
{
	this.p = leastSquaresForDiagonalIP(
		this.probabilityMatrix(), aimVector, this.nCandidates )

	pVector.fill(0)
	for (let i=0; i<this.p.length ; i++)	
	{	pVector.add(
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
		for (let i = 0; i< this.selection.length ; i++)
		{	this.p[i]+=(q[i])*adaptFactor
			
			// correct for overshoot
			if (this.p[i]<0)
			{	this.p[i]=0
			}
			pTotal+=this.p[i]
		}
		// message += "coeff : "+this.p
		console.log(message)
		if (null != htmlId)
		{	document.getElementById(htmlId).innerHTML = message
		}
		for (let i = 0; i< this.selection.length ; i++)
		{	
			// this.p[i]=this.p[i]/pTotal
		}
		// p = this.makePdistribution(p)
		adaptFactor *= stepFactor
	}
	return this.p
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

groups.prototype.pIndividuals = function(p,relative = false)
{	// p contains the probability of vectors in de selection
	var e = [];// expected number of groups that will have a place
	var a = []; // a[j]=e[j]/nGroups[j]; probability that a member of these groups will have a place
	for (let j = 0; j<this.groupSizes.length; j++)
	{
		e[j] = 0; // expected number of groups of groupSize[j]
		for (let i = 0; i< this.selection.length ; i++)
		{	e[j] += p[i]*this.combinations[this.selection[i]][j];
		}
		a[j]=e[j]/this.nGroups[j];
		if (relative ==true)
		{	a[j]=a[j]/this.pAim -1
		}
	}
	return a
}
groups.prototype.selectAllCombinations = function()
{	this.selection = Array.from(Array(this.combinations.length).keys())
}

groups.prototype.setCombinations = function()
{
	function recursivePart(nPlaces,nGroupsIn)
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
			{	return [[nPlaces/this.groupSizes[0]]]
			}
			// cannot fill all places without splitting a group
			return []
		}
		// else
		var result = []
			
		for (let i = 0; i<=max; i++)
		{	nLeft = nPlaces - i*this.groupSizes[nGroups.length]
			if (0<nLeft)
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
groups.prototype.calculationTeX = function(i, selec = this.selection)
{	
	if (selec.length <= i)
	{	return "There are not more than " + selec.length +" selected combinations."
			+ i + " is too big."
	}
	return this.calculationVectorToTeX ( this.combinations[selec[i]])
}

groups.prototype.calculationVectorToTeX = function(combi)
{	
	var result = "\\begin{pmatrix}"
	let total = 0
	for (j=0; j<this.groupSizes.length; j++)
	{	total +=  combi[j] * this.groupSizes[j]
		result += 0<j? "\\\\ " : ""
		result +=  "\\text{"+combi[j]+"x"+this.groupSizes[j] + "} &=& "
				+ combi[j] * this.groupSizes[j]
	}
	result += "\\\\ \\text{total} &=& " + total +"\\end{pmatrix} "
	
	return result 
}

groups.prototype.chancesTable = function(pVectors, nDigits = 3, pIndividuals = null)
{	
	if (null == pIndividuals)
	{	pIndividuals = this.pIndividuals(pVectors, false)
	}
	
	var result = "<table>"
	result +=    "<tr><th>Groups of size</th>"
		+    "<th>probability for persons<br>in this group </th><th>expected number <br> of groups of this size</th>"
	result +=    "<th>Expected number <br>of persons</th></tr>"

	var totalCandidates = 0
	var totalAdmitted = 0
	for (j=0; j<this.groupSizes.length; j++)
	{	let nCandidates = (this.nGroups[j]*this.groupSizes[j])
		let expectedGroups = (pIndividuals[j]*this.nGroups[j])
		let expectedAdmitted = expectedGroups*this.groupSizes[j]
		totalCandidates += nCandidates
		totalAdmitted += (expectedAdmitted)
		result += "<tr><td>" + this.nGroups[j]+" x " + this.groupSizes[j] 
			+ " person(s) = " + nCandidates
			+ "</td><td align='right'>"+pIndividuals[j].toFixed(8) +"</td>"
			+ "<td align='right'>" 
			+ this.nGroups[j] +	" x "+pIndividuals[j].toFixed(nDigits) 
			+ " = " 	+ expectedGroups.toFixed(nDigits)
			+ "</td><td align='right'>" +expectedGroups.toFixed(nDigits) + " x " 
			+ (this.groupSizes[j]) + " = "
			+ expectedAdmitted.toFixed(nDigits)
			+ "</td> </tr> \n"
	}
	result += "<tr><th align='right'> total: "+totalCandidates+"</th>" +
		  "<th align='right'>average P: "+this.pAim.toFixed(nDigits) + "</th><td></td>" +
		  "<th align='right'> total: "+totalAdmitted.toFixed(9) + "</th>" +
		"</table>"
	
	return result 
}
groups.prototype.resultCsv = function(selec = this.selection, p=this.p, separator = "\t", newline = "\r\n")
{
	result = ""
	// add headers
	result += "combination index-->"
	for (let i = 0; i<p.length; i++)
	{	result += separator + i
	}	
	// add Weights
	result += newline+"probability -->"
	for (let i = 0; i<p.length; i++)
	{	result += separator + p[i]
	}
	// add separator
	result += newline+"------ "
	for (let i = 0; i<p.length; i++)
	{	result += separator + "'-----'"
	}
	// add group header   
	result += newline+"groupsize " + separator + "selection -->"
	// add combinations
	for (j = 0; j<this.groupSizes.length; j++)
	{   result += newline + this.groupSizes[j]
	    for (let i = 0; i<p.length; i++)
	    {	result += separator + this.combinations[selec[i]][j]
	    }
	}
	// optional: add sums?
	return result
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
