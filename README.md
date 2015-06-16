# PowerSystemDataStructure.jl
Using a MATLAB power system data structure developed by a colleague, I crafted an equivalent Julia data structure.

## My dream data structure
I spent a year writing Visual Basic for Applications code. Though its syntax can be unwieldy, I found it ideal for working with Excel spreadsheets. The hierarchy of a power grid is similar to that of a spreadsheet. Ideally, I should be able to work with power system data using a VBA-esque syntax:
```vba
ActiveWorkbook.Sheets(1).Range("A1").Value
```
```julia
System.Lines[1].ThermalModel.mCp
```
This syntax is more human-readable than MATPOWER's caseformat.

In the last year I have written Julia code to perform instanton analysis and visualize the results. My code is messy. One method accepts a complete power system in mpc format, while another takes only a few pieces of data. This new power system data structure allows me to designate functions as system-level (accepting the top-level System object), line-level, etc.

Once data is loaded as a power system data instance, analysis can modify parameters and return an altered version of the same object. For example: transmission lines have no initial temperature data, but power flow may be used to assign a steady-state temperature to each line. Temperatures may be updated during subsequent dynamic analysis.
