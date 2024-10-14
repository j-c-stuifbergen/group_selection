/* e.g. 
nPlaces = 60;
groupSizes =[1,2,3,5]
nGroups = [60,8,2,1] 60 individuals, 8 duos, 2 trios, 1 group of 5
*/

function sumArray(q)
{	result = 0
	for (let i=0; i<q.length ; i++)
	{	result += q[i]
	}
	return result
}

/* make the total of p[] = 1
 * required: 0<=p[i] for each i
 * required: 0<sum of p[i]
 */
function makeDistribution(p)
{
	var pTotal =0;
	for (let i = 0; i< p.length ; i++)
	{	pTotal+=p[i];
	}	
	for (let i = 0; i< p.length ; i++)
	{	p[i]=p[i]/pTotal;
	}
	return p
}

function htmlFromArray(q)
{	return "[" + q + "]" 
}

function groups(nPlaces, groupSizes, nGroups) 
{
	this.nPlaces = nPlaces  // e.g. 80
	this.groupSizes = groupSizes // e.g. [1, 2, 3, 5] : groups can be 1, 2, 3 persons
	this.nGroups = nGroups  // e.g. [40, 8, 2, 1] for 40 individuals, 8 duos, 2 trios...
	this.nIndividuals = []
	this.combinations = [] // an array of combinations that fill all places
	this.selection = []
	this.p = [] // probabilities, belonging to combinations.
	this.pAim = this.averageP()
}


groups.prototype.individuals = function(groupSizes, nGroups)
{	result = new Array(groupSizes.length)
	for (let i=0; i<result.length ; i++)
	{	result[i] = groupSizes[i] * nGroups[i]
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
		  q[i]+=(this.pAim - a[j])*this.selection[i][j]*this.groupSizes[j];
		  // don't make chances negative
		  if (q[i] < 0 && 0 ===this.p[i])
		  {	q[i] = 0
		  }
		}
	}
	return q
}
groups.prototype.initialProbabilities = function()
{	
	this.p=[];
	for (let i = 0; i< this.selection.length ; i++)
	{	this.p[i]=1/this.selection.length;
	}
	
	return this.p
}
groups.prototype.improveProbabilities = function(nIterations, stepFactor = 1.001)
// at each iteration, adaptfactor will be multiplied by stepfactor
{
	if (null == this.p)
	{	this.p = this.initialProbabilities(this.selection)
	}
	var maxWeight = 0;
	for(i = 0; i<this.groupSizes.length; i++)
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

groups.prototype.averageP= function()
{
	// expected probability
	this.nCandidates = 0
	for (let i= 0; i<this.groupSizes.length; i++)
	{	this.nCandidates += this.groupSizes[i]*this.nGroups[i]
	}
	return this.nPlaces/this.nCandidates;
}
groups.prototype.pIndividuals = function(p,relative = false)
{
	var e = [];// expected number of groups that will have a place
	var a = []; // a[j]=e[j]/nGroups[j]; probability that a member of these groups will have a place
	for (let j = 0; j<this.selection[0].length; j++)
	{
		e[j] = 0; // expected number of groups of groupSize[j]
		for (let i = 0; i< this.selection.length ; i++)
		{	e[j] += p[i]*this.selection[i][j];
		}
		a[j]=e[j]/this.nGroups[j];
		if (relative ==true)
		{	a[j]=a[j]/this.pAim -1
		}
	}
	return a
}
groups.prototype.selectAllCombinations = function()
{	this.selection = this.combinations
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
groups.prototype.calculationString = function(i, combins = this.selection)
{	
	if (combins.length <= i)
	{	return "There are not more than " + combins.length +" selected combinations."
			+ i + " is too big."
	}
	var result = ""
	let total = 0
	let j = 0
	for (; j<this.groupSizes.length-1; j++)
	{	total +=  combins[i][j] * this.groupSizes[j]
		result += combins[i][j]+"x"+this.groupSizes[j]+" + "
	}
	total +=  combins[i][j] * this.groupSizes[j]
	result += combins[i][j]+"x"+this.groupSizes[j]+" = " + total
	
	return result 
}
groups.prototype.showCombinations = function(combins = this.selection, newLine = "<br>") // for console: "\r\n"
{
	let combinationString = ""
	let calcStr = ""
	for (let i = 0; i<combins.length; i++)
	{	combinationString += "["+combins[i]+"]"+newLine
		calcStr+= this.calculationString(i, combins)+newLine
 	}
	return {"combinations":combinationString, "calculations": calcStr}
}
