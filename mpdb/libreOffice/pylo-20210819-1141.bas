'mpdb: Material Property Data Base
'
'Revision Date: 2021.08.19
'
'SPDX-License-Identifier: BSD-2-Clause
'Copyright (c) 2021 Stuart Nolan. All rights reserved.
'
'Install:
'import this file into libreOffice
'
'create the directory ~/<path to libreOffice config>/<version>/user/Scripts/python
'and place python scripts (or links to them) in this directory
'e.g. On Ubuntu:
'lopyDir=${HOME}/.config/libreoffice/4/user/Scripts/python
'mpdbDir=${HOME}/local/src/github/mpdb/src/mpdb/
'mkdir ${lopyDir}
'cd ${lopyDir}
'ln -s ${mpdbDir}/libreOffice/pylo.py pylo.py
'ln -s ${mpdbDir}/units.py units.py
'ln -s ${mpdbDir}/periodicTable.py periodicTable.py
'ln -s ${mpdbDir}/costIndexes.py costIndexes.py

option explicit

Global masterScriptProvider as object

Function GetMasterScriptProvider() as object
	dim s as string
	dim factory as object
	s = "com.sun.star.script.provider.MasterScriptProviderFactory"
	if isNull(masterScriptProvider) then
		factory = createUnoService(s)
		masterScriptProvider = factory.createScriptProvider("")
	endif
	GetMasterScriptProvider = masterScriptProvider
End Function 

function GetScript(scriptName as string, scriptFile as string) as variant
	dim url as string
	dim msp as object
	url = "vnd.sun.star.script:" & scriptFile & "$" & scriptName & "?language=Python&location=user"
	msp = GetMasterScriptProvider()
	GetScript = msp.getScript(url)
end function

function cepci(year as variant, optional ciKey as variant, optional idxKey as variant) as variant
	dim script as object
	script = GetScript("cepci","costIndexes.py")
	if ismissing(idxKey) or isempty(idxKey) then
		idxKey = "year"
	end if
	if ismissing(ciKey) or isempty(ciKey) then
		ciKey = "ce plant cost index"
	end if
		
	cepci = script.invoke(Array(year,ciKey,idxKey), Array(), Array())
End Function

function cepci_keys() as variant
	dim script as object
	script = GetScript("cepci_keys","costIndexes.py")		
	cepci_keys = script.invoke(Array(), Array(), Array())
End Function

Function uc(value as double, inUnit as string, outUnit as string) as double
	dim script as object 
	script = GetScript("uc","units.py")
	uc = script.invoke(Array( value, inUnit, outUnit ), Array(), Array())
End Function

Function ucTemp(value as double, inUnit as string, outUnit as string) as double
	dim script as object
	script = GetScript("ucTemp","units.py")
	ucTemp = script.invoke(Array( value, inUnit, outUnit ), Array(), Array())
End Function

function fW(formula as string) as double
	dim script as object
	script = GetScript("fW","periodicTable.py")
	fW = script.invoke(Array(formula), Array(), Array())
End Function

function pySysVersion() as string
	dim script as object
	script = GetScript("pySysVersion","pylo.py")
	pySysVersion = script.invoke(Array(), Array(), Array())
End Function

function loArray(range as variant) as string
	dim script as object
	script = GetScript("loArray","pylo.py")
	loArray = script.invoke(Array(range), Array(), Array())
End Function

function interp(x as double, xVals as variant, yVals as variant, optional kind as variant) as variant
	dim script as object
	script = GetScript("interp","pylo.py")
	if ismissing(kind) or isempty(kind) then
		kind = "linear"
	end if
		
	interp = script.invoke(Array(x,xVals,yVals,kind), Array(), Array())
End Function

function interp2d(xvs as variant, yvs as variant, data as variant, xi as double, yi as double, optional meth as variant) as variant
	dim script as object
	script = GetScript("interp2d","pylo.py")
	if ismissing(meth) or isempty(meth) then
		meth = "linear"
	end if
		
	interp2d = script.invoke(Array(xvs, yvs, data, xi, yi, meth), Array(), Array())
End Function

function linInterp(x as double, xVals as variant, yVals as variant) as double
	'
	' TODO: input validtion:
	' 			numeric values only for x, xVals, and yVals, 
	'			blank values, 
	'			min(xVals) < x < max(xVals),
	'			extrapolation
	'
	Dim fAS : fAS = createUnoService( "com.sun.star.sheet.FunctionAccess" )
	dim x1i as integer
	dim x1 as double
	dim x2 as double
	dim y1 as double
	dim y2 as double

	x1i = fAS.CallFunction("MATCH", Array(x,xVals,1))
   	x1 = fAS.CallFunction("INDEX", Array(xVals,x1i))
   	y1 = fAS.CallFunction("INDEX", Array(yVals,x1i))
   	x2 = fAS.CallFunction("INDEX", Array(xVals,x1i+1))
   	y2 = fAS.CallFunction("INDEX", Array(yVals,x1i+1))
   	linInterp = y1 + (x - x1)*(y2 - y1)/(x2 - x1)
End function

