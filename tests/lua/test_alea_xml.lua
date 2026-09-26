local xml = [[
<alea version="1" length_units="cm" name="lua test">
  <materials>
    <material id="1" fraction_basis="atom">
      <nuclide zaid="1001" fraction="1" library="80c" />
    </material>
  </materials>
  <surfaces>
    <surface id="100000001" type="sphere" coeffs="0 0 0 2" />
  </surfaces>
  <cells>
    <cell id="1" name="inside" universe="0" material="1"
          density="1" density_units="g/cm3" region="-100000001">
      <importance particle="neutron" value="1" />
    </cell>
  </cells>
</alea>
]]

local system = alea.load_alea_string(xml)
assert(system:cell_count() == 1)
assert(system:surface_count() == 1)

local copy = system:export_alea_string()
assert(copy:find('<alea version="1"', 1, true))
assert(copy:find('library="80c"', 1, true))
assert(copy:find('particle="neutron"', 1, true))

local reloaded = alea.load_alea_string(copy)
assert(reloaded:cell_count() == 1)
print("test_alea_xml: OK")
