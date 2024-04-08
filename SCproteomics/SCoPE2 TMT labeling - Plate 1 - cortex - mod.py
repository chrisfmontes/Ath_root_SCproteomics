from opentrons import protocol_api

metadata = {
    'apiLevel': '2.14',
    'protocolName': 'SCoPE2 plate TMT 18-plex labeling - relabeling cortex',
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
    tiprack_2 = protocol.load_labware('opentrons_96_tiprack_20ul', 8)
    tiprack_3 = protocol.load_labware('opentrons_96_tiprack_20ul', 10)
    tiprack_4 = protocol.load_labware('opentrons_96_tiprack_20ul', 11)
    tiprack_1 = protocol.load_labware('opentrons_96_tiprack_20ul', 5)
    # Load single-chanel P20 on robot right arm
    p20 = protocol.load_instrument('p20_single_gen2', 'right', tip_racks=[tiprack_1, tiprack_2, tiprack_3, tiprack_4])

    #set temperature module at 10C and pause once it reaches temperature
    temperature_module.set_temperature(celsius=10)
    protocol.pause('Make sure TMT tubes and sample plate are in position and press "resume"')
    protocol.set_rail_lights(False)

    # We provide the position of wells to be labeled, one label at a time, from 5 (128C), to 18 (135N)
    myTMTdict = {
        "TMT5" :["C19","D7","D14","D24","E20","F13","G11","H16","H18","I3","I10","I18","J3","K5","L13","M3","N2","N13","N16","N17","O18","P10","P14"],
        "TMT6" :["F8","F16","G22","H7","H19","H21","J14","K17","L15","L23","M5","M12","N10","O14","O20","P3","P24"],
        "TMT7" :["E6","E8","E13","F17","G18","H11","I21","J8","J20","J23","J24","K1","K2","M6","M11","M22","N12","N18","N20","N23","O9","O17","O24"],
        "TMT8" :["E24","F3","H20","H23","I1","K4","K12","K14","K19","L6","L7","L8","L10","L21","M16","M20","M21","M23","N15","O1","O6","O12","P8","P12"],
        "TMT9" :["H2","H10","H12","H22","J15","J21","L20","M10","M15","N9","N24","O5","O8","P13","P19","P20"],
        "TMT10" :["F4","F7","G7","G17","I5","I12","I20","J9","K21","K24","N8","N11","N21","P4"],
        "TMT11" :["G19","H1","H17","H24","I9","I14","I17","I23","I24","J11","K6","K10","L12","M13","M18","N22","O7","O13","P5","P16"],
        "TMT12" :["D20","E19","E22","E23","F5","G9","G10","G12","G14","G23","G24","I4","J1","K16","M7","N1","O4","P7","P9","P22"],
        "TMT13" :["F11","F15","H3","H5","H8","H9","I13","J6","J19","J22","L5","L22","L24","O11","O19","O22","P15","P17","P21"],
        "TMT14" :["F10","F21","G13","G20","I6","I7","I15","I16","J4","J7","L1","L4","M2","M14","N4","O3","O10","O16","P6"],
        "TMT15" :["F18","F23","G4","H13","J2","J13","J16","K3","K7","K18","L9","L11","L18","L19","M17","N5","O2","O21","P1","P2"],
        "TMT16" :["G8","H6","J17","K11","K15","K20","L14","M1","M4","M9","N3","N14"],
        "TMT17" :["H14","I8","I11","J5","J18","K23","L2","L3","L17","M19","N7","N19","O15","P23"],
        "TMT18" :["H15","I2","I19","I22","J10","J12","K8","K9","K13","K22","L16","M8","M24","N6","O23","P11","P18"]
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