from opentrons import protocol_api

metadata = {
    'apiLevel': '2.14',
    'protocolName': 'SCoPE2 plate TMT 18-plex labeling',
    'author': 'Christian Montes - Walley Lab - Iowa State University',
    'description': 'Add TMTpro label to each well following provided well coordinates. Will prompt user to change label tube once current TMT label has been completely distributed. Magnetic module is not active and is only used to lift sample plate.'
    }

def run(protocol: protocol_api.ProtocolContext):
    protocol.set_rail_lights(True)
    # Load Magnetic Module GEN2 in deck slot 4 and sample plate on it.
    magnetic_module = protocol.load_module('magnetic module gen2', 4)
    sample_plate = magnetic_module.load_labware('thermofast_384_wellplate_40ul', label="SC samples")

    # Load Temperature Module GEN2 in deck slot 7 and 24-well aluminun block with 2mL tubes for TMT labels.
    temperature_module = protocol.load_module('temperature module gen2', 7)
    TMT_rack = temperature_module.load_labware('thermolowbindingtubes_24_aluminumblock_1500ul', label="TMT labels")

    #Load P20 tips
    tiprack_1 = protocol.load_labware('opentrons_96_tiprack_20ul', 8)
    tiprack_2 = protocol.load_labware('opentrons_96_tiprack_20ul', 10)
    tiprack_3 = protocol.load_labware('opentrons_96_tiprack_20ul', 11)
    tiprack_4 = protocol.load_labware('opentrons_96_tiprack_20ul', 5)
    # Load single-chanel P20 on robot right arm
    p20 = protocol.load_instrument('p20_single_gen2', 'right', tip_racks=[tiprack_1, tiprack_2, tiprack_3, tiprack_4])

    #set temperature module at 10C and pause once it reaches temperature
    temperature_module.set_temperature(celsius=10)
    protocol.pause('Make sure TMT tubes and sample plate are in position and press "resume"')
    protocol.set_rail_lights(False)

    # We provide the position of wells to be labeled, one label at a time, from 5 (128C), to 18 (135N)
    myTMTdict = {
        "TMT5": ["C11","C23","D2","D10","D15","E19","E21","E23","G11","G15","G17","H12","H15","I10","I19","J19","K3","K8","M1","M14","M23","N3","N10","N15","O7","O22"],
        "TMT6": ["A8","A15","A16","B7","C2","D9","D11","D16","D22","D24","E1","F17","F21","G9","G10","G14","H5","I1","I2","I22","J6","J10","J20","K20","L2","L3","L17","M17","N24","O15","P3","P11"],
        "TMT7": ["C17","D14","D19","D21","E3","E15","F1","F13","F24","G3","H21","I6","I11","I12","I16","I24","J14","K21","M7","M13","N9","O9","O16","P6","P16","P20"],
        "TMT8": ["B24","C4","C10","C16","D17","E16","F8","F20","G5","H1","H8","I4","I21","J3","K10","L1","L24","M4","O13","O19","P4","P5","P14","P15","P18"],
        "TMT9": ["A7","A13","B16","B23","C3","C18","D5","E8","E9","F7","F10","F12","G8","G16","H13","I17","I20","J5","J18","K2","K14","L10","N12","N20","N21","O10","P8","P10"],
        "TMT10": ["A19","B13","B14","C5","E2","F2","F5","F15","F16","F18","G7","G18","G22","G23","H18","I9","J22","J23","J24","K16","K17","K23","L6","L18","L19","N2","N13","N22","O2","O14"],
        "TMT11": ["A9","B17","C19","D23","E17","G1","G13","G21","H4","H6","H7","H24","I23","J15","K5","L8","L20","M6","M10","M18","N4","O6","P21","P22"],
        "TMT12": ["A14","B4","B18","B19","C7","C9","C20","C21","D20","E13","E14","H16","H19","I3","I18","J13","J21","K12","K19","M3","M9","M22","N5","O17"],
        "TMT13": ["A17","A20","B15","D4","E11","E18","F3","G20","I8","J8","J11","J12","J16","K4","K13","K22","L4","L15","M11","M20","N8","O1","O11","O24","P19"],
        "TMT14": ["A18","A23","B11","C24","D3","E7","E12","F6","F14","H11","I14","L5","L9","L21","M2","M12","M16","M19","M21","O3","O8","O12","O20","P7","P17"],
        "TMT15": ["A11","A21","B6","B12","C1","C8","C22","D8","D13","D18","E24","F22","G4","H22","I13","J2","J17","K1","K6","K11","N14","O21","P1","P23"],
        "TMT16": ["B1","B9","B22","C13","C14","C15","D7","D12","E20","F4","F9","H2","H9","I7","J1","J9","K18","L13","M5","M8","M15","M24","N1","N6","N17","N18","N19","O4","O5","O23","P9","P13"],
        "TMT17": ["A22","A24","B3","B5","B21","C12","D6","E5","E6","E10","G2","G6","G19","G24","H10","H14","H17","H23","I5","J4","K7","K9","K24","L12","L14","L22","L23","N11","P2","P24"],
        "TMT18": ["A10","A12","B2","B8","B10","B20","C6","D1","E4","E22","F11","F19","F23","G12","H3","H20","I15","J7","K15","L7","L11","L16","N7","N16","N23","O18","P12"]
        }
    for TMT in list(myTMTdict):
        p20.transfer(
        1,
        TMT_rack['D1'],
        [sample_plate.wells_by_name()[well_name] for well_name in myTMTdict[TMT]],
        new_tip='always',
        blow_out='True',
        blowout_location='destination well',
        mix_before=(1, 1))
        # And pause until next label is loaded into the aluminum block
        protocol.pause('Please replace TMT label tube on Temperature Module position D1, with next label')
        
    # p20.transfer(
        # 1,
        # TMT_rack['D1'],
        # [sample_plate.wells_by_name()[well_name] for well_name in ["A10","A16","B15","C6","C7","C11","D16","D17","E6","E12","E18","F23","G5","G16","H22","I5","I6","I10","I19","J14","K17","M16","M21","N5","N14","N16","O19","P5"]],
        # new_tip='once',
        # blow_out='True',
        # blowout_location='destination well',
        # mix_before=(1, 1))

    #deactivate temperature module at the end of protocol
    protocol.comment('Protocol finished, please remove plate and set it at room temperature for 2 hours')
    temperature_module.deactivate()