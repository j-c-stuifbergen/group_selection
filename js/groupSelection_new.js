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
function groups(nPlaces, groupSizes, nGroups) 
{
	this.nPlaces = nPlaces  // e.g. 80
	this.groupSizes = groupSizes // e.g. [1, 2, 3, 5] : groups can be 1, 2, 3 persons
	this.nGroups = nGroups  // e.g. [40, 8, 2, 1] for 40 individuals, 8 duos, 2 trios...
	this.nIndividuals = this.individuals(groupSizes, nGroups)
	
	// these two calculations should give the same result
	this.pAim = sumArray(this.nIndividuals)/this.nPlaces
	this.pAim = this.averageP()
}

groups.prototype.individuals = function(groupSizes, nGroups)
{	result = new Array(groupSizes.length)
	for (let i=0; i<result.length ; i++)
	{	result[i] = groupSizes[i] * nGroups[i]
	}
	return result
}

groups.prototype.derivativeIndividuals = function(p)
{	var q = []; // used to adapt p, store deviations from desired probabilities
	// calculate probabilities per individual 
	let a = this.pIndividuals(p)
	// calculate difference from desired probability
	for (let i = 0; i< this.combinations.length ; i++)
	{	q[i]=0;
		for (let j = 0; j<this.combinations[0].length; j++)
		{
		  // give extra weight if a deviation affects more persons
		  q[i]+=(this.pAim - a[j])*this.combinations[i][j];
		}
	}
	return q
}
groups.prototype.derivativeGroups = function(p)
{	var q = []; // used to adapt p, store deviations from desired probabilities
	// calculate probabilities per individual 
	let a = this.pGroups(p)
	// calculate difference from desired probability
	for (let i = 0; i< this.combinations.length ; i++)
	{	q[i]=0;
		for (let j = 0; j<this.combinations[0].length; j++)
		{
		  // give extra weight if a deviation affects more persons
		  q[i]+=(this.pAim - a[j])*this.combinations[i][j]*this.groupSizes[j];
		}
	}
	return q
}
groups.prototype.setProbabilities = function(nIterations,notNegP)
{
	if (null == notNegP)
	{	notNegP = 0.0001/this.combinations.length
	}

	var maxWeight = 0;
	for(i = 0; i<this.groupSizes.length; i++)
	{	maxWeight = Math.max(maxWeight,this.groupSizes[i]*this.nGroups[i])
	}

	var adaptFactor = 0.5/maxWeight/this.combinations.length
	var stepFactor = 1.001 // at each iteration, adaptfactor will be multiplied by this factor
	console.log("adaptFactor = "+adaptFactor)
		
	let p=[];
	// initialise probabilities 
	for (let i = 0; i< this.combinations.length ; i++)
	{	p[i]=1/this.combinations.length;
	}
	for (let k = 0; k< nIterations ; k++)
	{
		q= this.derivativeGroups(p)
		q2 = 0
		for(let i=0; i<q.length; i++)
		{	q2 += q[i]*q[i];
		}
		console.log("iteration "+k+", the norm of the derivative is "+Math.sqrt(q2))
		var pTotal = 0
		for (let i = 0; i< this.combinations.length ; i++)
		{	p[i]+=(q[i])*adaptFactor
			
			// correct for overshoot
			if (p[i]<0)
			{	p[i]=notNegP
			}
			pTotal+=p[i]
		}
		for (let i = 0; i< this.combinations.length ; i++)
		{	
			// p[i]=p[i]/pTotal
		}
		// p = this.makePdistribution(p)
		adaptFactor *= stepFactor
	}
	this.p = p
	return p
}

/* make the total of p[] = 1
 * required: 0<=p[i] for each i
 * required: 0<sum of p[i]
 */
groups.prototype.makePdistribution = function(p)
{
	var pTotal =0;
	for (let i = 0; i< this.combinations.length ; i++)
	{	pTotal+=p[i];
	}	
	for (let i = 0; i< this.combinations.length ; i++)
	{	p[i]=p[i]/pTotal;
	}
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
groups.prototype.pGroups = function(p,relative = false)
{
	var e = [];// expected number of groups that will have a place
	var a = []; // a[j]=e[j]/nGroups[j]; probability that a member of these groups will have a place
	for (let j = 0; j<this.combinations[0].length; j++)
	{
		e[j] = 0; // expected number of groups of groupSize[j]
		for (let i = 0; i< this.combinations.length ; i++)
		{	e[j] += p[i]*this.combinations[i][j];
		}
		a[j]=e[j]/this.nGroups[j];
		if (relative ==true)
		{	a[j]=a[j]/this.pAim -1
		}
	}
	return a
}
groups.prototype.pIndividuals = function(p,relative = false)
{
	var e = [];// expected number of persons that will have a place (sorted per groupSize)
	var a = []; // a[j]=e[j]/nGroups[j]; probability that a member of these groups will have a place
	for (let j = 0; j<this.combinations[0].length; j++)
	{
		e[j] = 0; // expected number of groups of groupSize[j]
		for (let i = 0; i< this.combinations.length ; i++)
		{	e[j] += p[i]*this.combinations[i][j];
		}
		a[j]=e[j]/(this.nGroups[j]*this.groupSizes[j]);
		if (relative ==true)
		{	a[j]=a[j]/this.pAim -1
		}
	}
	return a
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
groups.prototype.showCombinations = function(newLine = "<br>") // for console: "\r\n"
{
	let combinationString = ""
	let calculationString = ""
	for (let i = 0; i<this.combinations.length; i++)
	{	combinationString += "["+this.combinations[i]+"]"+newLine
		let total = 0
		let j=0
		for (j = 0; j<this.combinations[0].length-1; j++)
		{	calculationString += this.combinations[i][j]+"x"+this.groupSizes[j]+" + "
			total +=  this.combinations[i][j] * this.groupSizes[j]
		}
		total +=  this.combinations[i][j] * this.groupSizes[j]
		calculationString += this.combinations[i][j]+"x"+this.groupSizes[j]+" = " + total+newLine
 	}
	return {"combinations":combinationString, "calculations": calculationString}
}
