def main():
    import argparse

    from package.core import some_function  # Import your core functionality

    parser = argparse.ArgumentParser(description="My Python Package CLI")
    parser.add_argument("input", type=str, help="Input for the function")
    parser.add_argument("--option", type=str, help="Optional argument for the function")

    args = parser.parse_args()

    # Call the core function with the provided arguments
    result = some_function(args.input, args.option)
    print(result)


if __name__ == "__main__":
    main()
