import re
import csv
import logging

PFAM_FILE = r'Z:\users_files\Verena Burger\4_datasets\current\UniProt\pkinfam_n.csv'
UNIPROT_FILE = r"Z:\users_files\Verena Burger\4_datasets\current\UniProt\UniProt_human_reference_proteome.txt"
LOGGING_FILE = r'parseFlatFile.log'
logging.basicConfig(format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p', filename=LOGGING_FILE)


###
### definitions of functions
###
def pt(s):                                                              # where are these functions used?
	return "".join(["\t"]*s)

def oneInMany(source_range, target_range):
	for i in source_range:
		if i in target_range:
			return True

	return False

def expandY(seq, minStart = 0, maxEnd = 1000000):
	idx = seq.find("Y", minStart, maxEnd) + 1
	if idx == 0:
		return [seq]*2
	else:
		return [[seq[:idx], seq[:idx-1] + "pY"], expandY(seq[idx:], minStart = 0, maxEnd = min(maxEnd - idx, len(seq[idx:])))]

def collapseY(seqs):
	#print seqs
	if isinstance(seqs[1][0], str):
		if isinstance(seqs[1], str):
			return seqs
		return [a+b for a,b in zip(seqs[0], seqs[1])]
	else:
		ret = collapseY(seqs[1])
		return [a+b for a,b in zip([val for val in seqs[0] for _ in range(len(ret))], ret*len(ret))]

def getCombinatoricalpY(seq, minStart = 0, maxEnd = 10000000):
	return collapseY(expandY(seq, minStart = minStart, maxEnd = maxEnd))

def getSimplepY(seq, minStart = 0, maxEnd = 0):
	ls = minStart
	# special case for peptides which are at the beginning...
	if seq[ls-1:ls] == "Y" and seq[ls:ls+1] != "Y":
		ls -= 1
	if seq[ls:ls+1] == "Y":
		return [seq[:ls] + "pY" + seq[ls+1:]]
	else:
		return [seq]

def getYVersion(seq):
	ls = len(seq)/2
	return seq[:ls] + "Y" + seq[ls+1:]

def translateS(idx, ipos):
	ret = [i for i,x in enumerate(ipos) if x == idx]
	if len(ret) > 1:
		logging.info("ERROR translate, more than one return", idx, ipos, ret)
		exit(1)
	if len(ret) == 0:
		return None
	return ret[0]+1


def translate(idx, ipos):
	k = float(1)
	inc = 1
	K = float(10)
	start = int(len(ipos)*(k/K))-1
	while True:
		if isinstance(ipos[start], int):
			if ipos[start] < idx:
				k += inc
				start = int(len(ipos)*(k/K))-1
			else:
				break
		else:
			k += inc
		if k == K:
			start = len(ipos)-1
			break

	for i in range(start, -1, -1):
		if ipos[i] == idx:
			return i+1
	return None



#fopen = open("supplementary RTKs.txt", 'r')
#fopen = open("reviewed_rtk.txt", 'r')
#fopen = open("unreviewed_rtk.txt", 'r')
fopen = open(UNIPROT_FILE, 'r')

entry = {}
protein_list = {}
lastde = ""
type = ""
key = ""
ftid = 0
for line in fopen:
	line = line.rstrip()                                    # right strip: schneidet Leerzeichen rechts weg

	if line.startswith("ID"):
		logging.info("ID new id here")
		id_r = re.compile("ID\W*(\w*)\W*(\w*)\W*(\w*)")     # creates a regex named "id_r", notice the 3 groups!
		m = id_r.match(line)                                # matches (greps) the regex to the current ID line
		if m:
			entry["ID"] = {}                                # ? creating an empty object which is filled in the subsequent lines?
			entry["ID"]["NAME"] = m.group(1)                # backreferencing: extracting genename (probably^^)
			entry["ID"]["STATUS"] = m.group(2)              # extracting i.e. "reviewed"
		else:
			logging.info("ERROR ID", line)                  # no match

	if line.startswith("AC"):                               # I think those are the UniProt IDs of associated Proteins
		ac_r = re.compile("AC\W*(.*)")                      # making a regex
		m = ac_r.match(line)                                # matching the regex to the current line
		if m:
			if "AC" in entry:                               # schon wieder: wie kann etwas wo drin sein, das vorher als "leer" definiert wurde?
				te = m.group(1).replace(" ", "").split(";")[1:-1]  # concatenates the Uniprot IDs, separated by ","
				for e in te:
					entry["AC"].append(e)
			else:
				entry["AC"] = m.group(1).replace(" ", "").split(";")[0:-1]   # ? If there's just 1 UniProt ID?
		else:
			logging.info("ERROR AC", line)                  # no match

	if line.startswith("DE"):
		de_r = re.compile("DE\W*(.*)")                      # regex
		m = de_r.match(line)
		if m:
			if not("DE" in entry):
				entry["DE"] = {}
			de = m.group(1).split(": ")
			if len(de) == 2:
				entry["DE"][de[0]] = de[1]
				lastde = de[0]
			elif len(de) == 1:
				entry["DE"][lastde] = entry["DE"][lastde] + " " + de[0]
			else:
				logging.info("ERROR DE to much things", de)
			
		else:
			logging.info("ERROR DE", line)

	if line.startswith("GN"):
		gn_r = re.compile("GN\W*(.*)")
		m = gn_r.match(line)
		if m:
			gn = m.group(1).split(";")[0].split("=")
			if not("GN" in entry):
				entry["GN"] = {}
			if len(gn) == 3:
				gnt = gn[1].split(" ")
				entry["GN"][gn[0]] = " ".join(gnt[0:-1])[0:-1]
				entry["GN"][gnt[-1]] = gn[2][0:-1]
			elif len(gn) == 2:
				gnt = gn[1].split(" ")
				#entry["GN"][gn[0]] = gn[1][0:-1]
				entry["GN"][gn[0]] = gnt[0]
			else:
				logging.info("ERROR GN unknown amount", gn)
		else:
			logging.info("ERROR GN", line)

	if line.startswith("PE"):
		pe_r = re.compile("PE\W*(.*)")
		m = pe_r.match(line)
		if m:
			pe = m.group(1).split(":")
			entry["PE"] = pe[0]
		else:
			logging.info("ERROR PE", line)

	if line.startswith("FT"):
		ft_r1 = re.compile("FT\W+(\w+)\W+([0-9]+)\W+([0-9]+)\W*(.*)")
		ft_r2 = re.compile("FT\W*(.*)")
		m1 = ft_r1.match(line)
		m2 = ft_r2.match(line)
		if not("FT" in entry):
			entry["FT"] = {}
		if m1:
			if type != "" and key != "":
				if type in entry["FT"] and not("REPLACE" in entry["FT"][type][key]):
					repl_r = re.compile("([A-Z]+)\W*->\W*([A-Z]+)(.*)")
					replm_r = re.compile("Missing.*")
					isoform_r = re.compile(".*\(in\W*isoform\W*(.+)\).*")      # Wieso sucht er hier (in FT) nach den Isoformen?
					m = repl_r.match(entry["FT"][type][key]["ODESC"])
					mm = replm_r.match(entry["FT"][type][key]["ODESC"])
					iso = isoform_r.match(entry["FT"][type][key]["ODESC"])
					#if entry["ID"]["NAME"] == "FGFR2_HUMAN":
					#	logging.info "O",entry["FT"][type][key]["ODESC"] 
					if m:
						entry["FT"][type][key]["REPLACE"] = [m.group(1).replace(" ", ""), m.group(2).replace(" ", "")]
						if entry["FT"][type][key]["REPLACE"][1] == "Y":
							logging.info("THIS IS INTERESTING", type, key, entry["FT"][type][key])
							if not("YSPOT" in entry):
								entry["YSPOT"] = []
							entry["YSPOT"].append((type, key))
						elif not(type.startswith("VAR_SEQ")) and entry["FT"][type][key]["REPLACE"][1].find("Y") != -1:
							logging.info("THIS MIGHT ALSO BE INTERSTING", type, key, entry["FT"][type][key])

						entry["FT"][type][key]["DESC"] = m.group(3)
					elif mm:
						entry["FT"][type][key]["REPLACE"] = ["delete", ""]
						entry["FT"][type][key]["DELETE"] = ""
						entry["FT"][type][key]["DESC"] = entry["FT"][type][key]["ODESC"]
					else:
						entry["FT"][type][key]["DESC"] = entry["FT"][type][key]["ODESC"]
					if iso and (m or mm):
						isoformGroup = iso.group(1)
						isoforms = [c.split(",") for c in isoformGroup.replace("isoform", "").split("and")]
						isoforms = [isoform.strip() for sublist in isoforms for isoform in sublist]
						if not("ISOFORM" in entry):
							entry["ISOFORM"] = {}
						for isoform in isoforms:
							if not(isoform in entry["ISOFORM"]):
								entry["ISOFORM"][isoform] = []
							if mm:
								itype = "del"
							else:
								itype = "rep"
							entry["ISOFORM"][isoform].append((type, key, itype))
			ftid += 1
			type = m1.group(1) + "_" + str(ftid)
			start = m1.group(2)
			end = m1.group(3)
			desc = m1.group(4)
			if not(type in entry["FT"]):
				entry["FT"][type] = {}
			key = str(start)+"-"+str(end)
			if key in entry["FT"][type]:
				logging.info("ERROR FT already got this", type, key)
				exit(-1)
			entry["FT"][type][key] = {}
			entry["FT"][type][key]["ODESC"] = desc
		elif m2:
			if type in entry["FT"]:
				#key = str(start)+"-"+str(end)
				if key in entry["FT"][type]:
					entry["FT"][type][key]["ODESC"] = entry["FT"][type][key]["ODESC"] + m2.group(1)
				else:
					logging.info("ERROR FT no key", key)
			else:
				logging.info("ERROR FT no type", line, type)
		else:
			logging.info("ERROR FT", line)

	if line.startswith("SQ"):

		if type != "" and key != "":
			if "FT" in entry:
				entry["FT"][type][key]["DESC"] = entry["FT"][type][key]["ODESC"]
			type = ""
			key = ""

		sq_r = re.compile("SQ\W*SEQUENCE\W*([0-9]*) AA;\W*([0-9]*) MW;\W*(\w*) (\w*)")
		m = sq_r.match(line)
		if m:
			entry["SQ"] = {}
			#entry["SQ"]["AA"] = m.group(1)
			#entry["SQ"]["MW"] = m.group(2)
			#entry["SQ"][m.group(4)] = m.group(3)
			entry["SQ"]["SEQUENCE"] = ""
		else:
			logging.info("ERROR SQ", line)
		for line in fopen:
			line = line.strip()
			if line == "//":
				break
			entry["SQ"]["SEQUENCE"] = entry["SQ"]["SEQUENCE"] + line.replace(" ", "")


	if line == "//":
		if not(entry["ID"]["NAME"] in protein_list):        # makes a list of all entries in the reference proteome file, i.e.: RWDD1_HUMAN     Protein list ist doch bisher leer?!
			if "FT" in entry:
				entry["TAGS"] = {}
				for type in entry["FT"]:
					for key in entry["FT"][type]:
						if not(key in entry["TAGS"]):
							entry["TAGS"][key] = []
						entry["TAGS"][key].append((type,key))

			if "ISOFORM" in entry:# and entry["ID"]["NAME"] == "VGFR3_HUMAN":
				sequence = entry["SQ"]["SEQUENCE"]
				for isoform in entry["ISOFORM"]:
					mod_collection = []
					for mod in entry["ISOFORM"][isoform]:
						type = mod[0]
						key = mod[1]
						istart = int(key.split("-")[0]) - 1
						iend = int(key.split("-")[1])
						itype = mod[2]


						mod_collection.append((istart, iend, type, key, itype))
					mod_collection = sorted(mod_collection, key=lambda mod: -mod[0])

					isequence = sequence
					ipos = list(range(1, len(sequence)+2))
					#print isoform
					#print sequence
					#print len(sequence)
					for mod in mod_collection:
						istart = mod[0]
						iend = mod[1]
						type = mod[2]
						key = mod[3]
						itype = mod[4]
						if itype == "del":
							ipos = ipos[:istart] + ipos[iend:]
							#isequence = isequence[:istart] + "".join(["-"]*(iend-istart)) + isequence[iend:]
							isequence = isequence[:istart] + isequence[iend:]
							#logging.info istart, iend
							#logging.info ipos
							#logging.info isequence, len(isequence)
						elif itype == "rep":
							if not("REPLACE" in entry["FT"][type][key]):
								logging.info(entry["FT"][type][key])
							replace = entry["FT"][type][key]["REPLACE"]
							#logging.info replace, istart, iend
							if replace[0] == isequence[istart:iend]:
								ipos = ipos[:istart] + ["*"]*len(replace[1]) + ipos[iend:]
								isequence = isequence[:istart] + replace[1] + isequence[iend:]
							else:
								logging.info("ERROR replace does not match", replace[0], isequence[istart:iend])
								exit(-1)
							#logging.info ipos
							#logging.info isequence, len(isequence)
						else:
							logging.info("ERROR unkown itype", mod)
					
					entry["SQ"]["ISOFORM " + isoform] = (isequence, ipos, len(isequence))
				#logging.info entry["ID"]
				#logging.info entry["SQ"].keys()
				#logging.info entry["ISOFORM"]
			entry["SQ"]["SEQUENCE"] = (entry["SQ"]["SEQUENCE"], list(range(1, len(entry["SQ"]["SEQUENCE"]) + 2)), len(entry["SQ"]["SEQUENCE"]))   # Achtung: hier wird der Sequence name zusammengesetzt!!!
			protein_list[entry["ID"]["NAME"]] = entry                                           # liste der sequence names
		else:
			logging.info("ERROR protein list", entry["ID"]["NAME"])
		entry = {}
fopen.close()

kinase_group = dict()           # ein dictionary ist ein mapping file
with open(PFAM_FILE, 'r') as csvfile:
        kinases = csv.reader(csvfile, delimiter=',', quotechar='"')
        for kinase in kinases:
                if kinase[0] != "group" and kinase[3] != "":
                        if not(kinase[3] in kinase_group):
                                kinase_group[kinase[3]] = [kinase[0], kinase[1], kinase[2], kinase[4]]    # fills the dictionary with the kinase names (just human?) from the pkinfam file

for gene in list(protein_list.keys()):
	logging.info("GENE", gene)
	entry = protein_list[gene]
	sequence = entry["SQ"]["SEQUENCE"][0]
	#logging.info entry["FT"]
	has = 0
	if gene not in kinase_group:
		logging.info("ERROR, kinase", gene, "not in kinase_group")
		continue
	for type in list(entry["FT"].keys()):
		for key in list(entry["FT"][type].keys()):
			if type.startswith("DOMAIN"):
				logging.info(type, key, entry["FT"][type][key])
				desc = entry["FT"][type][key]["ODESC"]
				if "protein kinase" in desc or "Protein kinase" in desc or "PI3K/PI4K" in desc or "Histidine kinase" in desc:
					ystart = int(key.split("-")[0])
					yend = int(key.split("-")[1])
					egene = gene.split("_")[0]
					sgene = entry["GN"]["Name"]
					print("_".join(("".join((">", sgene)), kinase_group[gene][0], egene, key, kinase_group[gene][1], kinase_group[gene][2])))
					print(sequence[ystart-1:yend])
					has = has + 1
	logging.info(gene, "has", has)

