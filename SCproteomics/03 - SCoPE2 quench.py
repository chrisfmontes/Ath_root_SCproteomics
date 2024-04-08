from opentrons import protocol_api

metadata = {
    'apiLevel': '2.14',
    'protocolName': 'SCoPE2 TMT label quenching',
    'author': 'Christian Montes - Walley Lab - Iowa State University',
    'description': 'Add Hydroxylamine to TMT reacion and incubate to quench.'
    }

def run(protocol: protocol_api.ProtocolContext):
    # Turn OT-2 light on (because they are probably off)
    protocol.set_rail_lights(True)
  
    # Load Temperature Module GEN2 in deck slot 4 and 384-well sample plate on it.
    temperature_module = protocol.load_module('temperature module gen2', 4)
    sample_plate = temperature_module.load_labware('thermofast_384_wellplate_40ul', label="SC samples")
    
    # Load tube rack in slot 7 with 1.5 mL Eppendorf tube on it
    tube_rack = protocol.load_labware('opentrons_24_tuberack_eppendorf_1.5ml_safelock_snapcap',7)

    #Load P20 tips
    tiprack_2 = protocol.load_labware('opentrons_96_tiprack_20ul', 8)
    tiprack_3 = protocol.load_labware('opentrons_96_tiprack_20ul', 10)
    tiprack_4 = protocol.load_labware('opentrons_96_tiprack_20ul', 11)
    tiprack_1 = protocol.load_labware('opentrons_96_tiprack_20ul', 5)
    
    # Load single-chanel P20 on robot right arm
    p20 = protocol.load_instrument('p20_single_gen2', 'right', tip_racks=[tiprack_1, tiprack_2, tiprack_3, tiprack_4])

    # promt user to put sample plate and hydroxylamine tube
    protocol.pause('Make sure Hydroxylamine solution tube and sample plate are in position and press "resume"')
    temperature_module.set_temperature(celsius=10) #cool-down plate to prevent sample evaporation
    
    p20.transfer(
        1,
        tube_rack['D3'],
        sample_plate.wells(),
        new_tip='always',
        blow_out='True',
        blowout_location='destination well',
        mix_before=(1, 1)
        )
    
    # Incubate reaction
    protocol.pause('Please seal plate, then vortex, spin down and incubate 30 minutes at room temperature. Press "resume" to finish protocol.')
    temperature_module.deactivate()