#!/usr/bin/python3
#Synopsis:
#Here, all functions needed to execute the tasks of the program were collected. Different functions found here are called by the script 1_run_this_script.py depending on what task is selected. By user input.

#Importing what I need.
import os,shutil
from Bio import Entrez
#Should give a user email to Entrez, but would reveal identity:
Entrez.email = "Bxxxxxx-2020@ed.ac.uk"










#Function 1:
	#A function that asks for the input from the user to determine factors needed for making a blast database and running the blast.

#Function that downloads the fasta files:
def get_info(filename="123",db="protein",rmod="text",rtype="fasta"):
	from Bio import Entrez
	import os
	
	#This provides a dictionary containing all the Entrez databases
	databases = Entrez.einfo()
	databases_dict = Entrez.read(databases)	
	databases.close()

#Queries the user what database they would like to use and what search terms they want to employ
	db = input("Please select one of the following databases listed here and press enter.\n")
	print(databases_dict["DbList"],"\n\n")
	user_search_term = input("Thank you! \n Please provide search terms you would like to search the database for.\n")
	user_search = user_search_term.replace(" ","+")
	#os.system("clear")
	print ("Now searching the Entrez",db,"with the search term(s)", user_search)
#Error trap to ensure that the database inputted by the user is valid
	if db in databases_dict["DbList"]:
		print ("Your input database was valid")
	else:
		print ("An invalid database was entered. Returning to the interface.")
		go = input ("")
		os.system("./1_run_this_script.py")





#Searching Entrez db
	s_result = Entrez.esearch(db=db, term=user_search)
	record = Entrez.read(s_result)
	print (record["IdList"])

#Error trap to check if there were hits: 
	check = not bool(record["IdList"])
	if check == True:
		print ("Your search did not yield any results! Please press enter to return to the user interface and try again.")
		go = input("") 
		os.system("./1_run_this_script.py")

	print ("Got past error trap")




#This section gets the query key and webenv, which helps cope with bigger jobs 

	post = Entrez.read(Entrez.epost(db, id=",".join(record["IdList"])))
	webenv = post["WebEnv"]
	q_key = post["QueryKey"]
	#print (webenv)
	#print (q_key)

#Asking the user for parameters to run fasta download:
	filetype = input ("""Please enter what type of file you would like to download. This need to be a compatible abbreviation listed in the following table: \n
https://www.ncbi.nlm.nih.gov/books/NBK25499/table/chapter4.T._valid_values_of__e_and/?report=objectonly
\n\n The current default is the fasta file, but it is highly reccomended that the table comprised in this link is consulted when searching databases that are not the database protein or nuccore, as different parameters will have to be given successful.\n""")
	filemode = input ("""Please enter what format the file should be in. This needs to be a compatible abbreviation listed in the following table: \n
https://www.ncbi.nlm.nih.gov/books/NBK25499/table/chapter4.T._valid_values_of__retmode_and/?report=objectonly
\n\n The current default is xml, but it is highly reccomended that the table comprised in this link is consulted when searching databases that are not the database protein or nuccore, as different parameters will have to be given for the download to be successful.\n""")
	#Ask the user to specify a file name
	filename = input("Please provide a output file name.\n")
	#All the IDs obtained through querying the Entrez database specified by the user
	id_list=",".join(record["IdList"])

	return filename,db,rtype,rmod,id_list






#Function 2:
	#Function that downloads all the files from the Entrez databse using the parameters specified by the user in function 1 
def download(input_list):
	#If the filename does not exist, the function will cycle through all IDs and run an efetch on them
	if not os.path.isfile(input_list[0]):
		for ID in input_list[4:]:
			user_file=Entrez.efetch(db=input_list[1], id=ID, rettype=input_list[2], retmode=input_list[3])
			out_file=open(input_list[0],"w")
			out_file.write(user_file.read())
			out_file.close()
			user_file.close()
			print("Saved.")
	else:
		print("This file name already exists. Press enter to return to the interface.")
		next = input ("")
		os.system("./1_run_this_script.py")
	return input_list





#Running the functions 1 and 2:
#download_input = get_info()
#download(download_input)




#Function 3:
	#Making a blast database
def makedb (input_list):
	print ("generating database for the file:",input_list[0],",\nPlease wait...")

	#This generates a string that allows a system call to run makeblastdb. This was the easiest method to generate a blast database and all python methods I could find ended up using bash anyhow.
	#Asks what type of database you would like to generate
	db_type = input ("""Please provide the database type. In the current version, only the options nucl (nucleotide sequences database) and prot (protein sequences database) are supported. 
\n Please type in nucl to create a nucleotide database using your output file or prot to create a protein database using your output file.\n""")
	if db_type == "nucl" or db_type == "prot":
		bash_call = "makeblastdb -in "+input_list[0]+" -input_type  "+input_list[2]+" -dbtype "+db_type+" -parse_seqids -out  "+input_list[0]+".databasefile"
		
		os.system (str(bash_call))
		return input_list
	else:
		print ("Neither nucl nor prot were provided. Press enter to exit to the interface.\n")
		go =input ("")
		os.system("./1_run_this_script.py")
		return input_list









#Function 4:
	#Running a blast on the file created by the user
def blasting (input_list):
	from Bio.Blast import NCBIWWW
	from Bio import SeqIO
#Asking the user for what parameters are to be used in the blast.
	fasta_file = open(input_list[0]).read()
	blast_type = input("""Please specify what type of blast you would like to do. Currently supported versions are:
\nblastp : Used to compare protein sequences
\nblastn : Used to compare nucleotide sequences
\n\nKeep in mind that you will have needed to generate the corresponding database. \n
It is reccomended that the input file is a fasta file.\n """)
	out_form = input("""Please provide an output format for the BLAST:\n
HTML\n
Text\n
ASN.1\n
XML\n
""")	

#Error trapping for invalid inputs for the output format:
	if out_form == "HTML" or out_form == "Text" or out_form == "ASN.1" or out_form == "XML": 
		print (out_form,"Was selected")
	else:
		print ("An invalid input for selecting the output format was given. Please press enter to return to the interface.")
		go = input ("")
		os.system("./1_run_this_script.py")

#Error trapping for invalid input for the type of blast to be done:

#Blastp using user input	
	if blast_type == "blastp":
		database_type = input_list[0]+".databasefile"
		#print (database_type)
		print ("Doing blastp. Please wait...")
		from Bio.Blast.Applications import NcbiblastpCommandline
		blastp = NcbiblastpCommandline(query=input_list[0], db=database_type, outfmt=7, out=input_list[0]+"_blast")
		blastp()
		print ("A classic output format of 7 was selected for this blast to enable the processing of the data.")
	
#Blastn using user input
	elif blast_type == "blastn":
		from Bio.Blast.Applications import NcbiblastnCommandline
		database_type = input_list[0]+".databasefile"
		print (database_type)
		print ("Doing blastn. Please wait...")
		out_name = input_list[0]+"_blast"
		blastn = NcbiblastnCommandline(query=input_list[0], db=database_type, outfmt=7, out=out_name)
		print ("A classic output format of 7 was selected for this blast to enable the processing of the data.")
		blastn()
	#if (blast_type != "blastp") or (blast_type != "blastn"):
	#	print ("An invalid input for selecting the blast type was given. Please press enter to return to the interface.")
	#	go = input ("")
		#os.system("./1_run_this_script.py")
	print ("Blast is completed! Press enter to return to the interface.")
	go = input ("")	
	return input_list
















#Function 5:
	#Running a blast on whatever file the user desires. Much more breakable but allows for the user to enter any file containing proteins or sequences and blast them
def dbblasting (filename,inputtype):
	print ("generating database for the file:",filename,",\nPlease wait...")

	#This generates a string that allows a system call to run makeblastdb. This was the easiest method to generate a blast database and all python methods I could find ended up using bash anyhow.
	#Asks what type of database you would like to generate
	db_type = input ("""Please provide the database type. In the current version, only the options nucl (nucleotide sequences database) and prot (protein sequences database) are supported. 
\n Please type in nucl to create a nucleotide database using your output file or prot to create a protein database using your output file.\n""")
	if db_type == "nucl" or db_type == "prot":
		bash_call = "makeblastdb -in "+filename+" -input_type  "+inputtype+" -dbtype "+db_type+" -parse_seqids -out  "+filename+".databasefile"
		
		os.system (str(bash_call))
	else:
		print ("Neither nucl nor prot were provided. Press enter to exit to the interface.\n")
		go =input ("")
		#os.system("./1_run_this_script.py")
	from Bio.Blast import NCBIWWW
	from Bio import SeqIO
#Asking the user for what parameters are to be used in the blast.
	fasta_file = open(filename).read()
	blast_type = input("""Please specify what type of blast you would like to do. Currently supported versions are:
\nblastp : Used to compare protein sequences
\nblastn : Used to compare nucleotide sequences
\n\nKeep in mind that you will have needed to generate the corresponding database. \n
It is reccomended that the input file is a fasta file.\n """)
	out_form = input("""Please provide an output format for the BLAST:\n
HTML\n
Text\n
ASN.1\n
XML\n
""")	

#Error trapping for invalid inputs for the output format:
	if out_form == "HTML" or out_form == "Text" or out_form == "ASN.1" or out_form == "XML": 
		print (out_form,"Was selected")
	else:
		print ("An invalid input for selecting the output format was given. Please press enter to return to the interface.")
		go = input ("")
		os.system("./1_run_this_script.py")

#Error trapping for invalid input for the type of blast to be done:

#Blastp using user input	
	if blast_type == "blastp":
		database_type = filename+".databasefile"
		#print (database_type)
		print ("Doing blastp. Please wait...")
		from Bio.Blast.Applications import NcbiblastpCommandline
		blastp = NcbiblastpCommandline(query=filename, db=database_type, outfmt=7, out=filename+"_blast")
		blastp()
		print ("A classic output format of 7 was selected for this blast to enable the processing of the data.")
	
#Blastn using user input
	elif blast_type == "blastn":
		from Bio.Blast.Applications import NcbiblastnCommandline
		database_type = filename+".databasefile"
		print (database_type)
		print ("Doing blastn. Please wait...")
		out_name = filename+"_blast"
		blastn = NcbiblastnCommandline(query=filename, db=database_type, outfmt=7, out=out_name)
		print ("A classic output format of 7 was selected for this blast to enable the processing of the data.")
		blastn()
	#if (blast_type != "blastp") or (blast_type != "blastn"):
	#	print ("An invalid input for selecting the blast type was given. Please press enter to return to the interface.")
	#	go = input ("")
	#os.system("./1_run_this_script.py")
	print ("Blast is completed! Press enter to return to the interface.")
	go = input ("")	
	return filename,inputtype




