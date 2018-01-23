# VBT Changelog

### Version v1.0 :
	- VBT is now published at SBG Github!

### Version v0.7.7 :

	- Mendelian detailed logs table: variant category (SNP, INSERTION, DELETION) assign order changed from mother>father>child  to child>father>mother
	- Ref overlap mode becomes the default mode of VBT
	- Trim Variant from beginning is set to false by default in Mendelian tool

### Version v0.7.6 :

	- Quality column is added to the output trio of mendelian vcf
	- A new optional parameter is added that selects INFO tags from input file to be printed in the output trio file in Mendelian tool

### Version v0.7.5 :

	- VBT source is refactored (Duplicate code segment elimination, splitting long functions etc.) 
	- GPL License header is added to all files

### Version v0.7.4 :

	- Complex Skipped variant assigning is corrected (a lower boundary is introduced not to mark an included variant as complex skipped)
	- GA4GH output variant records with non BK/BD annotations are corrected

### Version v0.7.3 :
	
	- Sample test data is added under Example/
	- Readme file is modified
		- Detailed install guide is added
		- Sample execution command is added
	- General refactoring to increase code readability
