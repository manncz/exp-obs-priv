threshold_value = 20 
n = pd.read_csv(input_dat).shape[0]

describer = DataDescriber(category_threshold=threshold_value)

describer.describe_dataset_in_correlated_attribute_mode(dataset_file = input_dat, epsilon=syn_dp_ep, k=2)
describer.save_dataset_description_to_file(desc_file)

generator = DataGenerator()
generator.generate_dataset_in_correlated_attribute_mode(n, desc_file)

generator.save_synthetic_data(synth_dat)

synth_data = pd.read_csv(synth_dat)
