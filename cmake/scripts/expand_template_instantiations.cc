
#include <algorithm>
#include <cctype>
#include <fstream>
#include <iostream>
#include <ranges>
#include <regex>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

// NOLINTBEGIN

using TemplateExpansionMap = std::unordered_map<std::string, std::vector<std::string>>;

/**
 * A list of the available template expansions. The key is expansion list keyword, and the
 * accompanying list is the available values. For example, PRISMS_PF_DIMS -> (1, 2, 3)
 */
TemplateExpansionMap &
get_template_expansion_lists()
{
  static TemplateExpansionMap template_expansion_lists;
  return template_expansion_lists;
}

/**
 * Delete all whitespace in the front of a string
 */
void
delete_preceding_whitespace(std::string &str)
{
  auto it = std::ranges::find_if_not(str,
                                     [](unsigned char ch)
                                     {
                                       return std::isspace(ch);
                                     });
  str.erase(str.begin(), it);
}

/**
 * Delete all newlines
 */
void
delete_newlines(std::string &str)
{
  auto result = std::ranges::remove_if(str,
                                       [](char ch)
                                       {
                                         return ch == '\n' || ch == '\r';
                                       });
  str.erase(result.begin(), result.end());
}

/**
 * Delete comments
 */
void
delete_comments(std::string &str)
{
  static const std::regex single_line_comment {R"(//[^\n]*)"};
  static const std::regex multi_line_comment {R"(/\*[\s\S]*?\*/)"};

  str = std::regex_replace(str, multi_line_comment, "");
  str = std::regex_replace(str, single_line_comment, "");
}

/**
 * Delete a given set of characters in the front of a string
 */
void
delete_prefix_if_present(std::string &str, const std::string &prefix)
{
  if (str.starts_with(prefix))
    {
      str.erase(0, prefix.size());
    }
}

/**
 * Read the whole file into a single string variable.
 */
std::string
read_and_combine_file(std::istream &in_file)
{
  std::string file;
  std::string line;

  while (std::getline(in_file, line))
    {
      file += line + '\n';
    }

  return file;
}

/**
 * Extract the substring of the beginning of the provided string to the delimiter.
 */
std::string
get_substring_with_delimiter(std::string &str, const std::string &delimiter)
{
  std::string result;

  while (!str.empty())
    {
      char current = str.front();

      // Check if current char is in delimiter and not escaped
      bool is_delim = delimiter.find(current) != std::string::npos;
      bool escaped  = !result.empty() && result.back() == '\\';

      if (is_delim && !escaped)
        {
          break;
        }

      result += current;
      str.erase(str.begin());
    }
  // Trim trailing spaces from result
  while (!result.empty() && result.back() == ' ')
    {
      result.pop_back();
    }

  return result;
}

/**
 * Split a string into a vector with some delimiter
 */
std::vector<std::string>
split_string_to_vector(const std::string &str, const std::string &delimiter)
{
  std::vector<std::string> parts;
  size_t                   start = 0;
  size_t                   pos   = 0;

  while ((pos = str.find(delimiter, start)) != std::string::npos)
    {
      parts.push_back(str.substr(start, pos - start));
      start = pos + delimiter.size();
    }

  // Add the last part (or whole string if delimiter not found)
  // Make sure to clean spaces too
  auto part = str.substr(start);
  delete_preceding_whitespace(part);
  if (!part.empty())
    {
      parts.push_back(part);
    }

  return parts;
}

/**
 * Whether a string starts with a prefix
 */
bool
has_prefix(const std::string &str, const std::string &prefix)
{
  return str.starts_with(prefix);
}

/**
 * Replace a pattern into a string
 */
std::string
replace(const std::string &str,
        const std::string &pattern,
        const std::string &replacement)
{
  // If blank pattern return the original str
  if (pattern.empty())
    {
      return str;
    }

  std::string result = str;
  std::size_t pos    = 0;

  while ((pos = result.find(pattern, pos)) != std::string::npos)
    {
      result.replace(pos, pattern.length(), replacement);
      pos += replacement.length();
    }

  return result;
}

/**
 * Check that the provided token is the whole word and not part of something else. In
 * other words, check that only disallowed characters are on either side of the string
 * (e.g. &, <, etc..).
 */
bool
is_real_token(std::string_view       text,
              std::string::size_type pos,
              std::string::size_type length)
{
  constexpr std::string_view token_chars = "abcdefghijklmnopqrstuvwxyz"
                                           "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
                                           "0123456789_";

  auto is_token_char = [&](char c)
  {
    return std::ranges::find(token_chars, c) != token_chars.end();
  };
  if (pos > 0 && is_token_char(text[pos - 1]))
    {
      return false;
    }

  if (pos + length < text.size() && is_token_char(text[pos + length]))
    {
      return false;
    }

  return true;
}

/**
 * Replace a pattern into a string. This is a slightly more specialized version of the one
 * before since we pad the replacement with spaces along with some checks about word
 * connectivity.
 */
std::string
replace_token(std::string_view str,
              std::string_view pattern,
              std::string_view replacement)
{
  std::string result(str);
  std::size_t pos = 0;

  while ((pos = result.find(pattern, pos)) != std::string::npos)
    {
      if (is_real_token(result, pos, pattern.size()))
        {
          result.replace(pos,
                         pattern.size(),
                         std::string(" ") + replacement.data() + " ");
          pos += replacement.size() + 2;
        }
      else
        {
          ++pos;
        }
    }

  return result;
}

/**
 * Perform substitutions
 */
void
do_substitutions(
  const std::string                                      &str,
  const std::vector<std::pair<std::string, std::string>> &substitution_list)
{
  // We have to do this recursively to get all possible combinations. If the list is
  // empty, dump the text. Other perform the substitutions one at a time
  if (substitution_list.size() == 0)
    {
      std::cout << str << "\n" << std::flush;
      return;
    }

  // Grab the key value pairs
  const auto        &substitution_pair = substitution_list.back();
  const std::string &value             = substitution_pair.first;
  const std::string &key               = substitution_pair.second;

  // Check that the entry exists
  if (!get_template_expansion_lists().contains(key))
    {
      std::cerr << "Could not find a pattern for `" << key << "`\n" << std::flush;
      std::exit(1);
    }

  if (substitution_list.size() == 1)
    {
      // Loop over the expansion values for this key and recursively call this function
      for (const auto &expansion_value : get_template_expansion_lists()[key])
        {
          std::cout << replace_token(str, value, expansion_value) << std::endl;
        }
    }
  else
    {
      // Create a truncated list of substitutions, removing the back element
      const std::vector<std::pair<std::string, std::string>> remaining_substitution_list(
        substitution_list.begin(),
        substitution_list.end() - 1);

      // Loop over the expansion values for this key and recursively call this function
      for (const auto &expansion_value : get_template_expansion_lists()[key])
        {
          // Substitute the text with the expansion values
          std::string new_text = replace_token(str, value, expansion_value);

          // Perform substitutions
          do_substitutions(new_text, remaining_substitution_list);
        }
    }
}

/**
 * Read the template.in file
 */
void
read_templates(const std::string &filename)
{
  std::ifstream in_stream(filename);

  // Throw an error if we can't read the file
  if (!in_stream)
    {
      std::cerr << "Template file " << filename << " cannot be read\n" << std::flush;
      std::exit(1);
    }

  // Read the entire file and combine into a single string
  std::string file = read_and_combine_file(in_stream);

  // Delete the comments
  delete_comments(file);

  // Delete all newlines
  delete_newlines(file);

  // Remove the preceding spaces
  delete_preceding_whitespace(file);

  // Check that the string isn't empty. If it is we may have done something wrong with the
  // delete steps.
  if (file.empty())
    {
      std::cerr << "Template file " << filename << " is empty\n" << std::flush;
      std::exit(1);
    }

  // Process each line
  while (!file.empty())
    {
      // Grab the expansion name
      const std::string expansion_name = get_substring_with_delimiter(file, " :");

      // Remove the preceding spaces
      delete_preceding_whitespace(file);

      // Remove the :=
      delete_prefix_if_present(file, ":=");

      // Remove the preceding spaces
      delete_preceding_whitespace(file);

      // Check that next entry is a {, otherwise, the format is bad
      if (file[0] != '{')
        {
          std::cerr << "Invalid entry " << expansion_name << "\n" << std::flush;
          std::exit(1);
        }

      // Remove the {
      delete_prefix_if_present(file, "{");

      // Remove the preceding spaces
      delete_preceding_whitespace(file);

      // Grab the expansion values
      const std::string expansion_values = get_substring_with_delimiter(file, "}");

      // Check that next entry is a }, otherwise, the format is bad
      if (file[0] != '}')
        {
          std::cerr << "Invalid entry " << expansion_name << "\n" << std::flush;
          std::exit(1);
        }

      // Remove the }
      delete_prefix_if_present(file, "}");

      // Remove the preceding spaces
      delete_preceding_whitespace(file);

      // Assign the entry to the TemplateExpansionMap
      if (get_template_expansion_lists().contains(expansion_name))
        {
          std::cerr << "Invalid entry " << expansion_name << " already in the map\n"
                    << std::flush;
          std::exit(1);
        }
      get_template_expansion_lists()[expansion_name] =
        split_string_to_vector(expansion_values, ";");
    }
}

/**
 * Process the templates
 */
void
process_templates()
{
  // Grab the file the user has provided
  std::string class_template_in = read_and_combine_file(std::cin);

  // Delete comments
  delete_comments(class_template_in);

  while (!class_template_in.empty())
    {
      // Remove the preceding spaces
      delete_preceding_whitespace(class_template_in);

      // Keep preprocessor definitions
      if (has_prefix(class_template_in, "#"))
        {
          std::cout << get_substring_with_delimiter(class_template_in, "\n") << '\n';
          delete_preceding_whitespace(class_template_in);
          continue;
        }

      // Check that we have a for loop
      if (!has_prefix(class_template_in, "for"))
        {
          std::cerr << "Invalid template instantiation list. Missing `for`\n"
                    << std::flush;
          std::exit(1);
        }

      // Delete the for bit
      delete_prefix_if_present(class_template_in, "for");

      // Remove the preceding spaces
      delete_preceding_whitespace(class_template_in);

      // Check that we have a for loop
      if (!has_prefix(class_template_in, "("))
        {
          std::cerr << "Invalid template instantiation list. Missing `(`\n" << std::flush;
          std::exit(1);
        }

      // Delete the ( bit
      delete_prefix_if_present(class_template_in, "(");

      // Remove the preceding spaces
      delete_preceding_whitespace(class_template_in);

      // Grab the list of substitutions we have to make. These should correspond to the
      // expansion keys from read_templates()
      const std::vector<std::string> expansion_keys =
        split_string_to_vector(get_substring_with_delimiter(class_template_in, ")"), ";");

      // Check that we have a for loop
      if (!has_prefix(class_template_in, ")"))
        {
          std::cerr << "Invalid template instantiation list. Missing `)`\n" << std::flush;
          std::exit(1);
        }

      // Delete the ( bit
      delete_prefix_if_present(class_template_in, ")");

      // Remove the preceding spaces
      delete_preceding_whitespace(class_template_in);

      // Get a vector of all combinations of substitutions we have make
      std::vector<std::pair<std::string, std::string>> substitution_combinations;
      for (const auto &expansion_key : expansion_keys)
        {
          // Split the expansion keys
          const std::vector<std::string> value_key =
            split_string_to_vector(expansion_key, " :");
          if (value_key.size() != 2)
            {
              std::cerr << "Invalid template instantiation list. Expansions should be "
                           "paired key : value. You have `"
                        << expansion_key << "`\n"
                        << std::flush;
              std::exit(1);
            }
          const std::vector<std::string> values =
            split_string_to_vector(value_key[0], ",");

          std::transform(values.begin(),
                         values.end(),
                         std::back_inserter(substitution_combinations),
                         [&value_key](const auto &value)
                         {
                           return std::make_pair(value, value_key[1]);
                         });
        }

      // Now that we have the for loop bit sorted, read the interior

      // Remove the preceding spaces
      delete_preceding_whitespace(class_template_in);

      // Check that we have a for loop
      if (!has_prefix(class_template_in, "{"))
        {
          std::cerr << "Invalid template instantiation list. Missing `{`\n" << std::flush;
          std::exit(1);
        }

      // Delete the ( bit
      delete_prefix_if_present(class_template_in, "{");

      // Remove the preceding spaces
      delete_preceding_whitespace(class_template_in);

      // Grab the explicit template instantiations we have to substite
      const std::string substitution_templates =
        get_substring_with_delimiter(class_template_in, "}");

      // Delete the ( bit
      delete_prefix_if_present(class_template_in, "}");

      // Remove the preceding spaces
      delete_preceding_whitespace(class_template_in);

      // Replace the \{ with {
      const std::string new_substitution_templates =
        replace(replace(substitution_templates, "\\{", "{"), "\\}", "}");

      // Perform the substitutions
      do_substitutions(new_substitution_templates, substitution_combinations);
    }
}

int
main(int argc, char **argv)
{
  // Throw an error for the wrong script input
  if (argc < 2)
    {
      std::cerr << "Usage: \n  expand_template_instantiations list_of_in_files < inst.in "
                   "> inst \n"
                << std::flush;
      std::exit(1);
    }

  // Read the files that contain our templates
  for (int i = 1; i < argc; i++)
    {
      read_templates(argv[i]);
    }

  // Write a brief message for the .inst file
  std::cout << "// This file is automatically generated by PRISMS-PF\n\n" << std::flush;

  // Process the template instantiations
  process_templates();
}

// NOLINTEND
