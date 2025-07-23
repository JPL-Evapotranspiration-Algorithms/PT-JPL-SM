def verify() -> bool:
    import pandas as pd
    import numpy as np
    from .ECOv002_calval_PTJPLSM_inputs import load_ECOv002_calval_PTJPLSM_inputs
    from .process_PTJPLSM_table import process_PTJPLSM_table
    import os

    # Load input and output tables
    input_df = load_ECOv002_calval_PTJPLSM_inputs()
    module_dir = os.path.dirname(os.path.abspath(__file__))
    output_file_path = os.path.join(module_dir, "ECOv002-cal-val-PT-JPL-SM-outputs.csv")
    output_df = pd.read_csv(output_file_path)

    # Run the model on the input table
    model_df = process_PTJPLSM_table(input_df)

    # Columns to compare (model outputs)
    output_columns = [
        'G', 'Rn_soil', 'LE_soil', 'Rn_canopy', 'PET',
        'LE_canopy', 'LE_interception', 'LE'
    ]

    # Compare each output column
    for col in output_columns:
        if col not in model_df or col not in output_df:
            print(f"Missing column: {col}")
            return False
        # Use numpy allclose for floating point comparison
        if not np.allclose(model_df[col].values, output_df[col].values, rtol=1e-5, atol=1e-8, equal_nan=True):
            print(f"Mismatch in column: {col}")
            return False
    return True
