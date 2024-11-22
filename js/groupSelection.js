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
	let a = this.pIndividuals()
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
{	console.log("=== find the best combinations by least squares ===")

	aimVector = new Array(this.groupSizes.length).fill(this.pAim)
	console.log("aimVector " +aimVector)
	this.setProbabilities()

	pIndividu = Array(this.groupSizes.length).fill(0)

	this.selection = []
	this.p = []		
	
	this.selectFirstCombination() // perhaps not necessary?
//	this.p = [1] //otherwise, pIndividu will remain fill(0) 
	while (true)
	{	 
		pIndividu.fill(0)
		for (let i=0; i<this.p.length ; i++)	
		{	pIndividu = pIndividu.add(
			this.probabilities[ this.selection[i]]. timesScalar(this.p[i]))
		}

		pVerschil = aimVector.subtract(pIndividu)
		penalty = this.innerProd(pVerschil, pVerschil)
		// penalty = this.penalty(pIndividu)
		console.log("weights for the combinations are " + this.p)
		console.log("penalty =" + penalty + " ----- pVerschil =" + pVerschil)
		
		if ((this.p.length >= this.groupSizes.length)
		|| ( penalty < epsilon))
		{ 	break
		}
		ipResult = this.selectByIP(pVerschil)
		if (epsilon < ipResult.product )
		{	this.selection.push(ipResult.index)
		}
		else
		{	break
		}
		this.p = leastSquaresForDiagonalIP(
			this.probabilityMatrix(), aimVector, this.nIndividuals )
	}
			
	return this.makeDistribution() && (penalty < epsilon)
}

groups.prototype.selectByIP = function(goalVector)
{		 
	bestIndex = 0; biggestIP = 0
	let i = 0;
	for (; i<this.combinations.length ; i++)
	{	ip = this.innerProd(goalVector, this.probabilities[i])
		if (biggestIP < ip)
		{	biggestIP = ip
			bestIndex = i
		}
	}
	return {product: biggestIP, index:bestIndex}
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

groups.prototype.pIndividuals = function()
{	// p contains the probability of vectors in de selection
	var e = [];// expected number of groups that will have a place
	this.indivProba = []; // a[j]=e[j]/nGroups[j]; probability that a member of these groups will have a place
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
	}
	result += "<tr><th align='right'> total: "+totalCandidates+"</th>" +
		  "<th align='right'>average P: "+this.pAim.toFixed(nDigits) + "</th><td></td>" +
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
