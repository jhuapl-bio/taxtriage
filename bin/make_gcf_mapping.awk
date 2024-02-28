#!/usr/bin/awk -f

# Check if the line starts with ">"
/^>/ {
    line = substr($0, 2);

    # Find the first two spaces to separate chromosome and assembly from description
    firstSpace = index(line, " ");
    secondSpace = index(substr(line, firstSpace + 1), " ") + firstSpace;

    # Extract chromosome, assembly, and description
    chromosome = substr(line, 1, firstSpace-1);
    assembly = substr(line, firstSpace+1, secondSpace-firstSpace-1);
    description = substr(line, secondSpace+1);

    print chromosome "\t" assembly "\t" description;
}
