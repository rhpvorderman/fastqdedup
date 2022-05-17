#ifndef uint8_t
#include <stdint.h>
#endif
#ifndef size_t
#include <stddef.h>
#endif

static int 
within_hamming_distance(
    const uint8_t *string1,
    size_t string1_length,
    const uint8_t *string2,
    size_t string2_length,
    int max_distance)
{
    if (string1_length != string2_length) {
        // Hamming is technically only valid for sequences with the same
        // length. 
        return 0;
    }
    
    for (size_t i=0; i < string1_length; i++) {
        if (string1[i] != string2[i]) {
            max_distance -= 1;
            if (max_distance < 0) {
                return 0;
            }
        }
    }
    return 1;
}

static int
within_edit_distance(
    const uint8_t *string1,
    size_t string1_length,
    const uint8_t *string2,
    size_t string2_length,
    int max_distance)
{
    // Quick check for length differences.
    ssize_t length_difference = string1_length - string2_length;
    if (length_difference < 0) {
        length_difference = -length_difference;
    }
    if (length_difference > max_distance) {
        return 0;
    }
    while (string1_length > 0 && string2_length > 0) {
        if (string1[0] != string2[0]) {
            max_distance -= 1;
            if (max_distance < 0) {
                return 0;
            }
            int result;
            // Insertion. Compare string1's next character with current character string2.
            result = within_edit_distance(string1 + 1, 
                                          string1_length -1,
                                          string2,
                                          string2_length -1,
                                          max_distance);
            if (result) return result;
            // Deletion
            result = within_edit_distance(string1 + 1, 
                                          string1_length - 1, 
                                          string2 + 1,
                                          string2_length -1,
                                          max_distance);
            if (result) return result;
            // Otherwise we have a substitution. Continue the loop.
        }
        string1 += 1;
        string1_length -= 1;
        string2 += 1;
        string2_length -= 1;
    }
    // Strings are compared across their lengths, but if one is longer
    // than the other, it means we have a number of insertions the size of the 
    // difference.
    length_difference = string1_length - string2_length;
    if (length_difference < 0) {
        length_difference = -length_difference;
    }
    if (length_difference > max_distance) {
        return 0;
    }
    return 1;
}