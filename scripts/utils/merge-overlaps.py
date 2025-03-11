import pandas as pd

# Sample data
data = [
    [5332582, 5337837],
    [5360126, 5374508],
    [5375521, 5381384],
    [5381242, 5381384],
    [5381501, 5389833]
]

# Load into a DataFrame
columns = ["Ref_Start", "Ref_End"]
df = pd.DataFrame(data, columns=columns)

def remove_overlaps(df):
    """
    Removes overlapping intervals or merges them.
    Args:
        df (pd.DataFrame): DataFrame with 'Ref_Start' and 'Ref_End' columns.
    Returns:
        pd.DataFrame: Filtered DataFrame with no overlapping intervals.
    """
    # Sort intervals by start position
    df = df.sort_values(by="Ref_Start")

    # List to store non-overlapping intervals
    non_overlapping = []
    current_start, current_end = df.iloc[0]

    for i in range(1, len(df)):
        next_start, next_end = df.iloc[i]

        if next_start <= current_end:
            # Overlap detected; merge intervals
            current_end = max(current_end, next_end)
        else:
            # No overlap; save the current interval and move to the next
            non_overlapping.append([current_start, current_end])
            current_start, current_end = next_start, next_end

    # Add the last interval
    non_overlapping.append([current_start, current_end])

    # Convert back to DataFrame
    return pd.DataFrame(non_overlapping, columns=["Ref_Start", "Ref_End"])

# Apply the function
filtered_df = remove_overlaps(df)

# Display results
print(filtered_df)
