```yaml
# YAML CONFIG FILE
reference: "{path}/{to}/ref.fa"
overlap_regions: ["{path}/{to}/{tandem_repeat_regions}.bed"]
variant_sets:
    - type: "trCON"
      overlap_region_type: ["all"]
      repeat_count_change_range: [10, 50]
      number: 10
    - type: "trEXP"
      overlap_region_type: ["all"]
      repeat_count_change_range: [10, 50]
      number: 10
```