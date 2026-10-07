//! Every command shown in the user docs must use flags the real tools
//! accept. 0.5.2–0.5.8 documented `probe-demux --r2-fastq`,
//! `--probe-barcode-wl` and `--config`, none of which ever existed.
//!
//! Scans fenced code blocks in README.md, MANUAL.md and INSTALL.md, joins
//! `\`-continued lines, splits pipelines on `|`, and for each
//! `sbxpress <subcommand> ...` or `scripts/sbxpress*-pipeline.sh ...`
//! checks every `--flag` against that tool's `--help`.
use std::collections::BTreeMap;
use std::path::Path;
use std::process::Command;

const DOCS: &[&str] = &[
    "README.md",
    "MANUAL.md",
    "INSTALL.md",
    "docs/SERVER_TESTING.md",
];

fn code_block_commands(text: &str) -> Vec<(usize, String)> {
    let mut out = Vec::new();
    let mut in_block = false;
    let mut current = String::new();
    let mut start = 0;
    for (i, line) in text.lines().enumerate() {
        if line.trim_start().starts_with("```") {
            in_block = !in_block;
            current.clear();
            continue;
        }
        if !in_block {
            continue;
        }
        let line = match line.find(" #") {
            Some(p) if !line[..p].contains('"') => &line[..p],
            _ => line,
        };
        if current.is_empty() {
            start = i + 1;
        }
        let trimmed = line.trim_end();
        if let Some(body) = trimmed.strip_suffix('\\') {
            current.push_str(body);
            current.push(' ');
        } else {
            current.push_str(trimmed);
            out.push((start, std::mem::take(&mut current)));
        }
    }
    out
}

fn help_text(tool: &[&str], cache: &mut BTreeMap<String, String>) -> String {
    let key = tool.join(" ");
    cache
        .entry(key)
        .or_insert_with(|| {
            let mut cmd = if tool[0] == "sbxpress" {
                let mut c = Command::new(env!("CARGO_BIN_EXE_sbxpress"));
                c.args(&tool[1..]);
                c
            } else {
                let mut c = Command::new("bash");
                c.arg(Path::new(env!("CARGO_MANIFEST_DIR")).join(tool[0]));
                c
            };
            let out = cmd.arg("--help").output().unwrap();
            String::from_utf8_lossy(&out.stdout).to_string() + &String::from_utf8_lossy(&out.stderr)
        })
        .clone()
}

fn accepts(help: &str, flag: &str) -> bool {
    help.split(|c: char| {
        c.is_whitespace() || c == ',' || c == '[' || c == ']' || c == '/' || c == '=' || c == '<'
    })
    .any(|word| word == flag)
}

#[test]
fn documented_commands_use_real_flags() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR"));
    let mut cache = BTreeMap::new();
    let mut problems = Vec::new();
    let mut checked = 0;
    for doc in DOCS {
        let Ok(text) = std::fs::read_to_string(root.join(doc)) else {
            continue;
        };
        for (line, command) in code_block_commands(&text) {
            for segment in command.split('|') {
                let words: Vec<&str> = segment.split_whitespace().collect();
                let Some(pos) = words.iter().position(|w| {
                    *w == "sbxpress"
                        || w.trim_matches('"')
                            .ends_with("scripts/sbxpress-pipeline.sh")
                        || w.trim_matches('"')
                            .ends_with("scripts/sbxpress-gex-pipeline.sh")
                }) else {
                    continue;
                };
                let tool: Vec<&str> = if words[pos] == "sbxpress" {
                    match words.get(pos + 1) {
                        // minimap forwards flags to minimap2; bare --version/--help are global.
                        Some(sub) if !sub.starts_with('-') && *sub != "minimap" => {
                            vec!["sbxpress", sub]
                        }
                        _ => continue,
                    }
                } else {
                    let w = words[pos].trim_matches('"');
                    vec![&w[w.find("scripts/").unwrap()..]]
                };
                let help = help_text(&tool, &mut cache);
                checked += 1;
                for word in &words[pos + 1..] {
                    let flag = word.split('=').next().unwrap();
                    if flag.starts_with("--") && flag.len() > 2 && !accepts(&help, flag) {
                        problems.push(format!(
                            "{doc}:{line}: `{}` does not accept {flag}",
                            tool.join(" ")
                        ));
                    }
                }
            }
        }
    }
    assert!(
        checked > 10,
        "found only {checked} documented commands; is the parser broken?"
    );
    assert!(
        problems.is_empty(),
        "documented flags that do not exist:\n{}",
        problems.join("\n")
    );
}
