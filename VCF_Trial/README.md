# VCF File Consensus

Logan Wallace
10.15.2024

### Purpose

To create a 'high-confidence' set of structural variant calls for a patient by identifying a consensus between each of three different callers. 

### About

Hi all, in this folder are VCF files from three different callers for the same patient. This subset will serve as the trial data for our first attempt at identifying a consensus between the three callers. We will spike in some edge cases so we can be confident we are calling variants as we'd like. 

### Structure

We are starting with a single patient and a subset of their whole genome, but in the future will apply this to a cohort of approximately 60 patients. 

You'll need to make a comparison of the three different VCF files for this patient, identify which are concordant according to the rules below and then write them out to a final 'high confidence' VCF file which we will use for visualization in the second half of the project. 

### Rules for Consensus Identification

A variant should be written out to the consensus set if;

1. The variant is called in at least two of the three callers.
2. The variant must of the same variety (e.g., DEL) between the callers.
3. The variants must overlap each other to be called concordant.

### The Program

1. Do your best to make the program flexible. We may choose to change arguments for things like overlap or we may add an additional caller in the future. Keep an eye out for things you can allow the user to pass as input (argparse?) instead of hardcoding, etc. 
2. Annotate your code well and do your best to make it as readable as possible.
3. Try and make it portable. It would be great if this ends up being a 
4. Think about how you are approaching the problem and ways you might be able to restrict the search space to make the program more rapid or parallelizable. 

### Else

Take your time. We're hoping your program will be our tool for future problems of the same sort!

As always, please reach out if you have any questions! We are leaving a lot up to you all because you're smart and capable but don't suffer too long before asking for clarification. 

Please feel free to add to this github as you write code or notes or anything else you'd like to add. Now's a great time to practice some git skills. 

Best of luck and thank you!
