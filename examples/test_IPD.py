from CoDAC.InterproDomain import InterProDomain
import CODAC

# The lines below are used to generate the 'Domain_referencer_IPR000980_Human.csv' file
I1 = CODAC.Interpro('IPR000980', 'Homo sapiens')
uniprotlist = CODAC.Interpro.fetch_uniprotid(I1)
CODAC.Uniprot.makeRefFile(uniprotlist, 'IPR000980')

interpro_domain_class = InterProDomain()

metadata = interpro_domain_class.appendRefFile('Domain_referencer_IPR000980_Human.csv', 'Domain_referencer_IPR000980_Human_Interpro.csv')
