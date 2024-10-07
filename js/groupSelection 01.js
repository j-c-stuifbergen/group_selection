/* e.g. 
nPlaces = 60;
groupSizes =[1,40,50]
nGroups = [60,1,1]
*/
function groups(nPlaces, groupSizes, nGroups) 
{
	this.nPlaces = nPlaces  // e.g. 80
	this.groupSizes = groupSizes // e.g. [1, 2, 3, 5] : groups can be 1, 2, 3 persons
	this.nGroups = nGroups  // e.g. [40, 8, 2, 1] for 40 individuals, 8 duos, 2 trios...
}

groups.prototype.setProbabilities = function(nIterations,notNegP)
{
	if (null == notNegP)
	{	notNegP = 0.01/this.combinations.length
	}

	var maxWeight = 0;
	for(i = 0; i<this.groupSizes.length; i++)
	{	maxWeight = Math.max(maxWeight,this.groupSizes[i]*this.nGroups[i])
	}

	var adaptFactor = 0.5/maxWeight
	var stepFactor = 1.00 // at each iteration, adaptfactor will be multiplied by this factor
	console.log("adaptFactor = "+adaptFactor)
	var pAim = this.pAim()
		
	let p=[];
	// initialise probabilities 
	for (let i = 0; i< this.combinations.length ; i++)
	{	p[i]=1/this.combinations.length;
	}
	var q = []; // used to adapt p, store deviations from desired probabilities
	for (let k = 0; k< nIterations ; k++)
	{
		// calculate probabilities per individual 
		let a = this.pIndividuals(p)
		console.log("iteration "+k+", individual chances is "+a)
		// calculate deviations from desired probability
		qTotal = 0
		for (let i = 0; i< this.combinations.length ; i++)
		{	q[i]=0;
			for (let j = 0; j<this.combinations[0].length; j++)
			{
			  // size of the group is not taken into account
		//	  q[i]+=(pAim - a[j])*combinations[i][j];
			  // give extra weight if a deviation affects more persons
			  q[i]+=(pAim - a[j])*this.combinations[i][j]*this.groupSizes[j];
			  //q[i]+=(1/ a[j])*this.combinations[i][j]*this.groupSizes[j];
			}
			qTotal+=q[i]
		}
		var pTotal = 0
		for (let i = 0; i< this.combinations.length ; i++)
		{	
			// p[i]+=(q[i]-qTotal/this.combinations.length)*adaptFactor
			p[i]+=(q[i]/this.combinations.length)*adaptFactor
			// p[i]+=p[i]*(q[i]/this.combinations.length)*adaptFactor
			// p[i]+=(q[i])*adaptFactor
			
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


groups.prototype.pAim= function()
{
	// expected probability
	this.nCandidates = 0
	for (let i= 0; i<this.groupSizes.length; i++)
	{	this.nCandidates += this.groupSizes[i]*this.nGroups[i]
	}
	return this.nPlaces/this.nCandidates;
}
groups.prototype.pIndividuals = function(p,relative)
{
	var e = [];// expected number of groups that will have a place
	var a = []; // a[j]=e[j]/nGroups[j]; probability that a member of these groups will have a place
	var eTotal = 0
	var pAim = this.pAim()
	for (let j = 0; j<this.combinations[0].length; j++)
	{
		e[j] = 0; // expected number of groups of groupSize[j]
		for (let i = 0; i< this.combinations.length ; i++)
		{	e[j] += p[i]*this.combinations[i][j];
		}
		eTotal += e[j]*this.groupSizes[j]
		a[j]=e[j]/this.nGroups[j];
		if (relative ==true)
		{	a[j]=a[j]/pAim -1
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

	console.log("combinations"+this.combinations)
}
