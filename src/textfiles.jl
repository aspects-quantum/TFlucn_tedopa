
# Function to write outputs to a file with the current date and time
function write_to_file(filename::String, title::String, output::String)
    # Get the current date and time
    current_time = now()

    # Open the file in append mode
    open(filename, "a") do file
        # Write the current date and time to the file
        write(file, "Date and Time: $(current_time)\n")
        # Write the output to the file
        write(file, "Output: $(title)\n")
        write(file, "$(output)\n")
        # Add a separator for readability
        write(file, "----------------------\n")
    end
end

function write_for_loop(filename::String, i::String, output::String)
# Get the current date and time
current_time = now()
    # Open the file in append mode
    open(filename, "a") do file
        if i == string(1)
        # Add a separator for readability
        write(file, "----------------------\n")
        # Write the current date and time to the file
        write(file, ":::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::\n")
        write(file, "                                     Date and Time: $(current_time)\n")
        write(file, ":::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::\n")
        # Write the output to the file
        write(file, "$(output);\n")
        else
            write(file, "$(output);\n")
        end
    end
end

