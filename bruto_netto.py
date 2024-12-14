# dit programma dient om een grafiek te maken dat het verband tussen
# bruto-inkomen en netto-inkomen weergeeft.

def procent(percentage):
	return percentage*0.01

class tarief:
    def __init__(self, naam):
    	self.schijven = []
    	self.naam = naam

    def voegSchijfToe(self, vanaf, percentage, offset = None):
    	if (None == offset):
            if (0==len(self.schijven)):
                offset = 0;
            else:
                index = len(self.schijven)-1
		offset = tarief.berekenInSchijf(vanaf, self.schijven[index])
                
        self.schijven.append({ "vanaf": vanaf,
                        "perc": percentage*procent(1),
                        "offset": offset
                        })

    @staticmethod
    def vergelijk(schijf1, schijf2):
       return schijf1[vanaf] < schijf2[vanaf]

    @staticmethod
    def berekenInSchijf(schijf, inkomen):
        return ((inkomen - (schijf["vanaf"]))
                        * schijf["perc"]
                        + schijf["offset"])

    def berekening(self, inkomen):
        schijfIndex = len(self.schijven)-1
        while (inkomen < self.schijven[schijfIndex].vanaf):
           schijfIndex -= 1
        return tarief.berekenInSchijf(inkomen, self.schijven[schijfIndex])
    
    @staticmethod
    def inkomstenBelasting():
    	result = tarief("inkomstenBelasting")
	result.voegSchijfToe(0,39.93,0)
	result.voegSchijfToe(73032,49.5)
	
	return result

    @staticmethod
    def algemeneHeffingsKorting():
    	result = tarief("AlgemeneHeffingskorting")
	result.voegSchijfToe(0,0,-3070)
	result.voegSchijfToe(22661,6.095,-3070)
	result.voegSchijfToe(73031,0,0)
	
	return result

    @staticmethod
    def arbeidsKorting():
    	result = tarief("arbeidsKorting")
	result.voegSchijfToe(0,-8.23,0)
	result.voegSchijfToe(10741,-29.86,884)
	result.voegSchijfToe(23200,-3.09,4605)
	result.voegSchijfToe(37691,+6.51,5052)
	result.voegSchijfToe(115295,0,0)

	return result

    @staticmethod
    def bijdrageZorgverzekeringsWet():
    	result = tarief("zorgverzekeringsWet")
	result.voegSchijfToe(0,5.43,0)
	result.voegSchijfToe(66956,0,None)
	
	return result

class brutoNetto:
    def __init__(self):
	self.tarieven = [] # van class tarief natuurlijk
    
    def voegTariefToe(self, tarief):
    	if ("tarief" == typeOf(tarief) ):
	    tarieven.append(self.tarief)
	else:
	    throwException("invoer is niet van class 'tarief'")

    def abscissaas(self):
    	# use a set in stead of an array, to avoid duplicate values
	abscis = {} 
	for tarief in self.tarieven:
	    for schijf in tarief.schijven:
		abscis.add(schijf.vanaf)
	# make a list so the results can be sorted
	result = list(abscis).sort()
	result.append(result[len(result)-1] * 1.2)

	return result

    def heffing(self, inkomen):
    	naam = "naam"
	bedragTag = "bedrag"
	nettoPosition = 1

	resultaat = [ {naam: 'bruto',
			bedragTag: inkomen
			}]
	totaal = 0;
	netto = bruto;

	for tarief in self.tarieven:
	    heffing = tarief.berekening(inkomen)
	    netto -= heffing
	    resultaat.append( {naam:tarief.naam,
  				bedragTag: heffing
			})

	resultaat.insert(nettoPosition, {naam:'netto',
  					bedragTag: netto}
			)
	return resultaat

    def csvVanTarieven(self):
	gegevens = self.heffing(self.abscissaas)
	resultaat = ""

	for geg in gegevens:
	    resultaat += geg["naam"] + csvDelimiter
	resultaat += csvEoL
	
	for i in range(len(gegevens[0]["bedrag"])):
		for geg in gegevens:
		    resultaat += geg["bedrag"][i] + csvDelimiter
		resultaat += csvEoL

	return resultaat

    def maakTarieven2023(self):
	self.tarieven.append(tarief.inkomstenBelasting())	
	self.tarieven.append(tarief.bijdrageZorgverzekeringsWet())	
	self.tarieven.append(tarief.algemeneHeffingsKorting())	
	self.tarieven.append(tarief.arbeidsKorting())	

mijnBN = brutoNetto()
mijnBN.maakTarieven2023()

print (mijnBN.csvVanTarieven())
