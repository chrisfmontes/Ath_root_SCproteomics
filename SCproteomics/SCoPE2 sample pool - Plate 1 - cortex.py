from opentrons import protocol_api

metadata = {
    'apiLevel': '2.14',
    'protocolName': 'SCoPE2 TMT sample pooling',
    'author': 'Christian Montes - Walley Lab - Iowa State University',
    'description': 'Pool TMT-labeled and quenched samples into autosampler glass inserts for MS runs.'
    }

def run(protocol: protocol_api.ProtocolContext):
    # Turn OT-2 light on (because they are probably off)
    protocol.set_rail_lights(True)
  
    # Load Temperature Module GEN2 in deck slot 4 and 384-well sample plate on it.
    temperature_module = protocol.load_module('temperature module gen2', 4)
    sample_plate = temperature_module.load_labware('thermofast_384_wellplate_40ul', label="SC samples")
    
    # Load tube rack in slot 7 with 1.5 mL Eppendorf tube on it
    tube_rack = protocol.load_labware('opentrons_24_tuberack_eppendorf_1.5ml_safelock_snapcap',7)

    # Load glass vial adapter plate into slot 1
    glassVials = protocol.load_labware('thermoautosamplerglassinsertvialonfoamadapter_96_wellplate_250ul',1)
    
    #Load P20 tips
    tiprack_2 = protocol.load_labware('opentrons_96_tiprack_20ul', 8)
    tiprack_3 = protocol.load_labware('opentrons_96_tiprack_20ul', 10)
    tiprack_4 = protocol.load_labware('opentrons_96_tiprack_20ul', 11)
    tiprack_1 = protocol.load_labware('opentrons_96_tiprack_20ul', 5)
    tiprack_5 = protocol.load_labware('opentrons_96_tiprack_20ul', 2)
    
    # Load single-chanel P20 on robot right arm
    p20 = protocol.load_instrument('p20_single_gen2', 'right', tip_racks=[tiprack_1, tiprack_2, tiprack_3, tiprack_4, tiprack_5])

    # promt user to put sample plate and Optima H2O tube
    protocol.pause('Make sure sample plate and glass vials are in position and press "resume"')
    temperature_module.set_temperature(celsius=12) #cool-down plate to prevent sample evaporation
     
    # Define wells per plex
    myTMTplex = {
        "PlexNumber1" :["A24","D13","E5","E7","L2","M22","P9"],
        "PlexNumber2" :["B6","C18","F7","I6","K12","K17","O23"],
        "PlexNumber3" :["B1","B9","E20","K20","L20","M17","N22"],
        "PlexNumber4" :["H10","M13","N3","N6","N20","O15","P15"],
        "PlexNumber5" :["A8","A19","C23","I12","L15","N12","O19"],
        "PlexNumber6" :["A17","I2","I7","J13","K2","P3","P4"],
        "PlexNumber7" :["A15","E12","G10","O10","O24","P2","P18"],
        "PlexNumber8" :["D12","D22","E10","I8","I19","N8","P13"],
        "PlexNumber9" :["C17","E4","F23","H7","I16","M1","P8"],
        "PlexNumber10" :["B16","B22","D8","D11","K7","M24","P12"],
        "PlexNumber11" :["D3","F3","I4","J10","M11","N11","O8"],
        "PlexNumber12" :["D19","F2","F17","J17","M14","O20","P7"],
        "PlexNumber13" :["B14","C13","G6","G15","K22","L7","P19"],
        "PlexNumber14" :["B24","C16","G2"], #removed G7, C5, F15, and L17
        "PlexNumber15" :["D21","G18","I5","I17","K19","L3","O3"],
        "PlexNumber16" :["B5","G7","F5","F21","K23","M18","O9"],
        "PlexNumber17" :["D5","E11","F11","I23","J3","M5","M21"],
        "PlexNumber18" :["G14","H8","K1","K3","K13","O6","O18"],
        "PlexNumber19" :["F15","C5","E2","F10","G8","H2","I14"],
        "PlexNumber20" :["B2","B10","M2","M6","M19","N5","P5"],
        "PlexNumber21" :["A9","F19","F20","J24","K18","N1","O16"],
        "PlexNumber22" :["C1","C3","F22","I10","K10","M23","P17"],
        "PlexNumber23" :["E19","F6","H17","I13","L16","N16","P20"],
        "PlexNumber24" :["E1","E9","I21","L4","L10","M3","P21"],
        "PlexNumber25" :["A12","A20","B17","H3","H14","H22","J8"],
        "PlexNumber26" :["C14","E16","H23","I20","J15","J20","J22"],
        "PlexNumber27" :["E21","E23","I15","L9","M15","M20","P14"],
        "PlexNumber28" :["B21","E22","G17","I9","J7","K8","L13"],
        "PlexNumber29" :["A10","C20","J12","K16","L14","L24","N10"],
        "PlexNumber30" :["C11","C19","D17","E17","H24","J18","L22"],
        "PlexNumber31" :["A14","C6","C8","E18","F14","K24","P16"],
        "PlexNumber32" :["B18","E13","F9","G12","G13","H13","O11"],
        "PlexNumber33" :["C21","D2","G21","I24","L8","M7","N9"],
        "PlexNumber34" :["D1","H18","K4","L23","N19","O22","P22"],
        "PlexNumber35" :["B3","D16","D20","I22","J5","K11","O17"],
        "PlexNumber36" :["B19","C22","G11","H9","L11","M16","O7"],
        "PlexNumber37" :["A13","A22","C15","D7","J14","L19","N4"],
        "PlexNumber38" :["A16","G1","G24","I18","M9","N7","O12"],
        "PlexNumber39" :["B20","L17","F18","G22","I1","J6","N14"],
        "PlexNumber40" :["B15","H11","H20","J9","J19","L12","O5"],
        "PlexNumber41" :["D4","F13","H4","J21","J23","O2","P23"],
        "PlexNumber42" :["B8","C10","D9","D24","H21","L18","M10"],
        "PlexNumber43" :["E6","G19","H5","J1","K9","K15","N13"],
        "PlexNumber44" :["B7","B11","C4","D15","F8","G3","P6"],
        "PlexNumber45" :["C2","F4","F16","G5","J2","J4","P11"],
        "PlexNumber46" :["G20","I11","L5","L6","N21","O21","P10"],
        "PlexNumber47" :["A18","E8","I3","K14","L1","O4","P1"],
        "PlexNumber48" :["D6","D23","G23","H15","H19","M4","N24"],
        "PlexNumber49" :["C12","E14","E15","H12","J16","N2","N18"],
        "PlexNumber50" :["C7","F24","G16","H1","N17","O1","O14"],
        "PlexNumber51" :["A21","A23","D14","E24","H6","K6","K21"],
        "PlexNumber52" :["B23","D18","E3","H16","L21","O13","P24"],
        "PlexNumber53" :["D10","F1","G4","G9","J11","M8","N15"],
        "PlexNumber54" :["A7","A11","C24","F12","K5","M12","N23"]
        }
        
    destCoord = ["A1","A2","A3","A4","A5","A6","A7","A8","A9","A10","A11","A12","B1","B2","B3","B4","B5","B6","B7","B8","B9","B10","B11","B12","C1","C2","C3","C4","C5","C6","C7","C8","C9","C10","C11","C12","D1","D2","D3","D4","D5","D6","D7","D8","D9","D10","D11","D12","E1","E2","E3","E4","E5","E6"]
    
    for PLEX, WELL in zip(myTMTplex, destCoord):
        p20.consolidate(
        5,
        [sample_plate.wells_by_name()[well_name] for well_name in myTMTplex[PLEX]],
        glassVials[WELL]
        )
    

    # Incubate reaction
    protocol.pause('Please seal plate, then vortex, spin down and incubate 30 minutes at room temperature. Press "resume" to finish protocol.')
    temperature_module.deactivate()