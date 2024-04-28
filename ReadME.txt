Please verify the SHA256 Hash: AA8862BB53498FA4B32E9DF2E737B59C02C35B0A26C27438E7F9B9585BD36FA5 before running the executable to ensure the file hasn't been tampered with
Do this by opening PowerShell window in dist folder and running the command: <(Get-FileHash Blowdown_Analysis.exe -Algorithm SHA256).Hash>  

RUN APPLICATION (after SHA-256 Checksum verification): Open dist folder and run Blowdown_Analysis.exe 
USAGE INSTRUCTIONS: https://www.youtube.com/watch?v=aYoRy9T3CfQ
SOURCE CODE: Run Blowdown_Analysis.py
NOTE: Blowdown_Analysis.exe is only supported on Windows. For Mac, clone the Repo and run the python source code.

BACKGROUND:

This program calculates pressure, temperature and flow rate of a real gas over time during a
blowdown process. A blowdown process is a depressurization event where gas escapes from a
pressurized volume through a metering device (typically a valve, orifice, nozzle or leak).

Two blowdown processes are analyzed: isothermal and adiabatic. The difference between these
processes lies in how temperature and heat transfer are managed during the depressurization
of the gas.


