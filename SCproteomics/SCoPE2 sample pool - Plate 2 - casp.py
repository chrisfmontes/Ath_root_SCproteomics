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
        "PlexNumber1" :["C2","D13","G1","G8","J7","O12","C4"],
        "PlexNumber2" :["A9","A11","B18","E8","F24","J9","K9"],
        "PlexNumber3" :["E4","E11","G19","K12","O8","O14","P4"],
        "PlexNumber4" :["A14","C5","D3","D8","H12","L1","M17"],
        "PlexNumber5" :["E24","F23","M2","N4","N12","O5","O7"],
        "PlexNumber6" :["C16","F12","G13","J19","K7","M22","N17"],
        "PlexNumber7" :["B17","P15","C15","D24","G11","J8","P10"],
        "PlexNumber8" :["B12","C14","D17","H19","I16","I23","M1"],
        "PlexNumber9" :["B20","H24","K8","K23","K24","O10","O11"],
        "PlexNumber10" :["B19","J10","N2","N8","N15","O9","P13"],
        "PlexNumber11" :["A18","B22","E23","G14","L20","L22","N14"],
        "PlexNumber12" :["C1","G21","J24","L10","L16","M4","N10"],
        "PlexNumber13" :["A15","C24","F13","H18","I10","J2","K13"],
        "PlexNumber14" :["B10","C3","E15","K1","M14","P3","P5"],
        "PlexNumber15" :["A13","A21","B4","L4","M23","N24","P12"],
        "PlexNumber16" :["B16","B24","E22","H5","M5","O22","P19"],
        "PlexNumber17" :["D14","F15","H9","H16","J4","J18","O18"],
        "PlexNumber18" :["F10","F21","H4","M16","N19","O2","P2"],
        "PlexNumber19" :["C23","F2","H3","J14","J20","L23","O19"],
        "PlexNumber20" :["A7","D4","D10","F9","G9","K17","N5"],
        "PlexNumber21" :["A8","C19","F16","K22","N20","O13","P24"],
        "PlexNumber22" :["B2","G2","G4","I2","L6","P9","P16"],
        "PlexNumber23" :["D12","F1","L24","N11","N22","O15","O21"],
        "PlexNumber24" :["A19","B1","D16","F11","K6","M18","O17"],
        "PlexNumber25" :["A10","A23","C7","F17","H15","J22","K18"],
        "PlexNumber26" :["A16","A22","E21","F14","J13","K11","N6"],
        "PlexNumber27" :["B14","C6","C13","E5","F3","I1","N9"],
        "PlexNumber28" :["C10","H10","I13","J5","M7","O23","P11"],
        "PlexNumber29" :["B5","C8","G16","H7","L18","P18","P20"],
        "PlexNumber30" :["D1","D5","D18","I9","I22","M3","O16"],
        "PlexNumber31" :["A24","B11","D9","D19","F8","G15","K4"],
        "PlexNumber32" :["B9","B23","C11","D6","F20","J23","L17"],
        "PlexNumber33" :["A17","B6","E2","G17","H11","J1","L7"],
        "PlexNumber34" :["H2","I12","N7","O3","O6","P8","P23"],
        "PlexNumber35" :["D11","E7","G7","H1","H6","K3","L15"],
        "PlexNumber36" :["E12","G6","I3","L2","M13","N1","N13"],
        "PlexNumber37" :["E10","E13","I21","K14","L11","M15","M20"],
        "PlexNumber38" :["B7","F6","G12","H21","J16","P1","P22"],
        "PlexNumber39" :["C9","F7","F19","I24","K5","N3","P7"],
        "PlexNumber40" :["D15","F4","F22","H14","J6","L21","N16"],
        "PlexNumber41" :["E14","G10","G23","K10","L8","O1","O4"],
        "PlexNumber42" :["A20","E3","G5","J15","J21","M12","N23"],
        "PlexNumber43" :["C22","E1","H23","I20","L5","L19","P14"],
        "PlexNumber44" :["E16","G24","H20","I17","J12","K16","K19"],
        "PlexNumber45" :["C12","D2","D7","E17","G3","H8","J11"],
        "PlexNumber46" :["A12","C20","D21","D22","E20","M10","N21"],
        "PlexNumber47" :["E9","G22","I5","I8","K15","N18","P21"],
        "PlexNumber48" :["B3","E19","G18","I11","L9","M6","M11"],
        "PlexNumber49" :["C21","I4","I15","K20","L12","M24","P17"],
        "PlexNumber50" :["B15","C17","C18","F18","I18","J17","M21"],
        "PlexNumber51" :["B8","B21","D20","E18","I14","L3","P6"],
        "PlexNumber52" :["G20","H22","I6","K2","L14","M8","M9"],
        "PlexNumber53" :["B13","E6","I19","K21","L13","O20","O24"],
        "PlexNumber54" :["D23","F5","H13","H17","I7","J3","M19"]
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