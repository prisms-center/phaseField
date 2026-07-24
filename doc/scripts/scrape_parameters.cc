#include <deal.II/base/parameter_handler.h>

#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>

#include <prismspf/user_inputs/user_input_parameters.h>

#include <filesystem>
#include <fstream>
#include <string>

struct ParameterNode
{
  std::string            section;
  std::string            name;
  std::string            default_value;
  std::string            pattern_description;
  std::list<std::string> aliases;
  std::string            documentation;
  std::string            warnings;
  std::string            example;
  bool                   deprecated = false;
};

class DocumentationParser
{
public:
  struct ParsedDocumentation
  {
    std::string description;
    std::string warnings;
    std::string example;
  };

  static ParsedDocumentation
  parse(const std::string &raw_doc)
  {
    constexpr std::string_view warning_tag     = "@warning";
    constexpr std::string_view warning_end_tag = "@end_warning";
    constexpr std::string_view example_tag     = "@example";
    constexpr std::string_view example_end_tag = "@end_example";

    ParsedDocumentation result;

    const std::size_t warning_pos = find_tag(raw_doc, std::string(warning_tag));
    const std::size_t example_pos = find_tag(raw_doc, std::string(example_tag));

    const std::size_t first_tag_pos = std::min(warning_pos, example_pos);
    result.description              = trim(raw_doc.substr(0, first_tag_pos));

    if (warning_pos != std::string::npos)
      {
        const std::size_t content_start = warning_pos + warning_tag.size();

        const std::size_t end_marker_pos =
          find_tag(raw_doc, std::string(warning_end_tag));

        const std::size_t warning_end =
          (end_marker_pos != std::string::npos) ? end_marker_pos : raw_doc.size();

        result.warnings =
          extract_section(raw_doc, warning_pos, content_start, warning_end);
      }

    if (example_pos != std::string::npos)
      {
        const std::size_t content_start = example_pos + example_tag.size();

        const std::size_t end_marker_pos =
          find_tag(raw_doc, std::string(example_end_tag));

        const std::size_t example_end =
          (end_marker_pos != std::string::npos) ? end_marker_pos : raw_doc.size();

        result.example =
          extract_section(raw_doc, example_pos, content_start, example_end);
      }

    return result;
  }

private:
  static std::size_t
  find_tag(const std::string &doc, const std::string &tag)
  {
    return doc.find(tag);
  }

  static std::string
  extract_section(const std::string &doc,
                  std::size_t        tag_start,
                  std::size_t        content_start,
                  std::size_t        section_end)
  {
    if (tag_start == std::string::npos || content_start == std::string::npos)
      {
        return "";
      }

    const std::string raw = doc.substr(content_start, section_end - content_start);
    return trim(raw);
  }

  static std::string
  trim(const std::string &str)
  {
    const auto start = str.find_first_not_of(" \t\n\r");
    if (start == std::string::npos)
      {
        return "";
      }

    const auto end = str.find_last_not_of(" \t\n\r");
    return str.substr(start, end - start + 1);
  }
};

class ParameterParser
{
public:
  explicit ParameterParser(const std::string &file)
  {
    // Check that the file exists
    std::ifstream _file(file);
    AssertThrow(_file.is_open(),
                dealii::ExcMessage("Could not open parameter file: " + file));
    _file.close();

    // Read the boost property tree
    boost::property_tree::ptree root;
    boost::property_tree::read_json(file, root);

    // Walk the property tree, collecting nodes
    flatten(root, "");

    // Resolve aliases from collected data
    resolve_aliases();
  }

  [[nodiscard]] const std::vector<ParameterNode> &
  get_parameters() const
  {
    return _parameters;
  }

private:
  void
  flatten(const boost::property_tree::ptree &node, const std::string &prefix)
  {
    for (const auto &[key, child] : node)
      {
        // Skip if the entry is empty
        if (key.empty())
          {
            continue;
          }

        const std::string full_path = prefix.empty() ? key : prefix + "." + key;

        if (is_parameter_entry(child))
          {
            ParameterNode param;

            param.section             = prefix.empty() ? "root" : prefix;
            param.name                = key;
            param.default_value       = child.get<std::string>("default_value", "");
            param.pattern_description = child.get<std::string>("pattern_description", "");
            param.deprecated          = child.get<bool>("deprecation_status", false);

            const auto raw_doc = child.get<std::string>("documentation", "");
            const DocumentationParser::ParsedDocumentation parsed =
              DocumentationParser::parse(raw_doc);

            param.documentation = parsed.description;
            param.warnings      = parsed.warnings;
            param.example       = parsed.example;

            _parameter_index[full_path] = _parameters.size();

            _parameters.push_back(std::move(param));
          }
        else if (is_alias_entry(child))
          {
            const std::string canonical_name = child.get<std::string>("alias", "");

            if (!canonical_name.empty())
              {
                const std::string canonical_path =
                  prefix.empty() ? canonical_name : prefix + "." + canonical_name;

                _raw_aliases[full_path] = {canonical_path, key};
              }
          }
        else if (!child.empty())
          {
            flatten(child, full_path);
          }
      }
  }

  [[nodiscard]] bool
  is_parameter_entry(const boost::property_tree::ptree &node) const
  {
    // deal.II parameter entries always have a value and default entry. Subsections don't,
    // hence the need for this function
    return node.find("value") != node.not_found() &&
           node.find("default_value") != node.not_found();
  }

  [[nodiscard]] bool
  is_alias_entry(const boost::property_tree::ptree &node) const
  {
    // deal.II alias entries always have a leaf with only a alias entry.
    return node.find("alias") != node.not_found();
  }

  void
  resolve_aliases()
  {
    for (const auto &[alias_path, alias_info] : _raw_aliases)
      {
        const auto &[canonical_path, alias_name] = alias_info;

        auto canonical_it = _parameter_index.find(canonical_path);

        if (canonical_it == _parameter_index.end())
          {
            auto redirect_it = _raw_aliases.find(canonical_path);
            if (redirect_it != _raw_aliases.end())
              {
                canonical_it = _parameter_index.find(std::get<0>(redirect_it->second));
              }
          }

        ParameterNode &target = _parameters[canonical_it->second];
        target.aliases.push_back(alias_name);
      }
  }

  std::vector<ParameterNode> _parameters;

  std::map<std::string, std::pair<std::string, std::string>> _raw_aliases;

  std::map<std::string, std::size_t> _parameter_index;
};

class HTMLGenerator
{
public:
  explicit HTMLGenerator(const std::vector<ParameterNode> &parameters)
    : _parameters(parameters)
  {}

  void
  write(const std::string &input_filepath, const std::string &output_filepath) const
  {
    // Read the input file
    std::ifstream in(input_filepath);
    AssertThrow(in.is_open(),
                dealii::ExcMessage("Could not open file: " + input_filepath));

    std::ostringstream buffer;
    buffer << in.rdbuf();
    std::string content = buffer.str();
    in.close();

    // Find the marker
    const std::string marker     = "@parameter_list";
    const std::size_t marker_pos = content.find(marker);
    AssertThrow(marker_pos != std::string::npos,
                dealii::ExcMessage("Could not find marker " + marker + " in file " +
                                   input_filepath));

    // Build the replacement content
    const std::string sections = build_sections();

    // Replace the marker with the generated sections
    content.replace(marker_pos, marker.size(), sections);

    // Write to output
    std::ofstream out(output_filepath);
    AssertThrow(out.is_open(),
                dealii::ExcMessage("Could not open file: " + output_filepath));
    out << content;
    out.close();
  }

private:
  [[nodiscard]] std::map<std::string, std::vector<const ParameterNode *>>
  group_by_section() const
  {
    std::map<std::string, std::vector<const ParameterNode *>> groups;

    for (const auto &param : _parameters)
      {
        groups[param.section].push_back(&param);
      }

    // Sort parameters within each section alphabetically
    for (auto &[section, params] : groups)
      {
        std::sort(params.begin(),
                  params.end(),
                  [](const ParameterNode *a, const ParameterNode *b)
                  {
                    return a->name < b->name;
                  });
      }

    return groups;
  }

  [[nodiscard]] std::string
  build_sections() const
  {
    const auto sections = group_by_section();

    std::ostringstream ss;

    // Builds root section first and without title
    auto root_it = sections.find("root");
    if (root_it != sections.end())
      {
        ss << build_table(root_it->second);
      }

    for (const auto &[section, params] : sections)
      {
        if (section == "root" || section == "Root")
          {
            continue;
          }
        ss << build_section(section, params);
      }

    return ss.str();
  }

  [[nodiscard]] std::string
  build_table(const std::vector<const ParameterNode *> &params) const
  {
    std::ostringstream ss;
    // Open the HTML table
    ss << "<table>\n"
       << "<tr>\n"
       << "<th>Parameter</th>\n"
       << "<th>Default</th>\n"
       << "<th>Pattern</th>\n"
       << "<th>Description</th>\n"
       << "<th>Aliases</th>\n"
       << "</tr>\n";

    for (const auto *param : params)
      {
        ss << build_row(*param);
      }

    ss << "</table>\n\n";

    return ss.str();
  }

  [[nodiscard]] std::string
  build_section(const std::string                        &section,
                const std::vector<const ParameterNode *> &params) const
  {
    std::ostringstream ss;

    // Doxygen section heading
    ss << "\\section " << to_anchor(section) << " " << section << "\n\n";

    ss << build_table(params);

    return ss.str();
  }

  [[nodiscard]] std::string
  build_row(const ParameterNode &param) const
  {
    std::ostringstream ss;

    // Deprecated rows get an inline style for the red highlight.
    const std::string row_style = param.deprecated
                                    ? " style=\"background-color: rgba(192,57,43,0.15);"
                                      " border-left: 4px solid #c0392b;\""
                                    : "";

    // Name the cell based on the canonical parameter name. If the parameter is
    // deprecated, add on a deprecated status.
    std::ostringstream name_cell;
    name_cell << "<span class=\"parameter-badge\" onclick=\"copy_parameter(this)\">"
              << param.name << "</span> ";
    if (param.deprecated)
      {
        name_cell << " <b>(deprecated)</b>";
      }

    // Now we move onto the description, which combines the documentation, warnings, and
    // examples
    std::ostringstream desc_cell;
    desc_cell << param.documentation;

    if (!param.warnings.empty())
      {
        desc_cell << "<br>@warning " << param.warnings << "\n";
      }

    if (!param.example.empty())
      {
        desc_cell << "<br>@note Here's an example of how it would look in the "
                     "parameters.prm file";
        desc_cell << "@code\n" << param.example << "\n@endcode";
      }

    // Now we move onto aliases. We format these as little badges since there is usually
    // of lot of them.
    std::ostringstream alias_cell;
    if (!param.aliases.empty())
      {
        for (const auto &alias : param.aliases)
          {
            alias_cell << "<span class=\"alias-badge\" onclick=\"copy_parameter(this)\">"
                       << alias << "</span> ";
          }
      }

    // Now we onto the patterns
    // TODO: Format patterns to be a little more readable

    ss << "<tr" << row_style << ">\n"
       << "<td>" << name_cell.str() << "</td>\n"
       << "<td><tt>" << param.default_value << "</tt></td>\n"
       << "<td>" << param.pattern_description << "</td>\n"
       << "<td>" << desc_cell.str() << "</td>\n"
       << "<td>" << alias_cell.str() << "</td>\n"
       << "</tr>\n";

    return ss.str();
  }

  static std::string
  to_anchor(const std::string &section)
  {
    std::string anchor = section;

    // replace spaces, dots, and colons with underscores, lowercase everything
    std::transform(anchor.begin(),
                   anchor.end(),
                   anchor.begin(),
                   [](const unsigned char c) -> char
                   {
                     if (c == ' ' || c == '.' || c == ':')
                       {
                         return '_';
                       }
                     return (char) std::tolower(c);
                   });

    // also prefix the section label with p_list_
    return "p_list_" + anchor;
  }

  const std::vector<ParameterNode> &_parameters;
};

int
main(int argc, char *argv[])
{
  // Throw an error for the wrong script input
  if (argc != 3)
    {
      std::cerr << "Usage: \n  scrape_parameters input_file output_file\n" << std::flush;
      std::exit(1);
    }

  // Check that the parameters file exists
  const std::string in_file = argv[1];
  AssertThrow(std::filesystem::exists(std::filesystem::path(in_file)),
              dealii::ExcMessage("Could not find file " + in_file));
  const std::string out_file = argv[2];

  // Create the deal.II parameter handler and the PRISMS-PF
  prismspf::UserInputParameters<3> user_inputs;
  dealii::ParameterHandler         parameter_handler;

  // Print the parameter handler object to json
  std::ofstream json_file("parameters.json");
  user_inputs.declare(parameter_handler, 1);
  parameter_handler.print_parameters(json_file, dealii::ParameterHandler::JSON);

  // Read back in the parameter into a more useful structure
  ParameterParser parser("parameters.json");
  const auto      parameters = parser.get_parameters();

  // Write to HTML
  HTMLGenerator writer(parameters);
  writer.write(in_file, out_file);

  return 0;
}
