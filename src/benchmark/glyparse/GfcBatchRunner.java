import java.io.BufferedReader;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.nio.charset.StandardCharsets;
import java.util.Base64;

import org.glycoinfo.GlycanFormatConverter.cli.ConverterPortal;
import org.glycoinfo.GlycanFormatConverter.cli.InputFormat;
import org.glycoinfo.GlycanFormatConverter.cli.OutputFormat;

/** Persistent-JVM adapter for the glyparse corpus benchmark. */
public final class GfcBatchRunner {
    private static final String PROTOCOL = "glyparse-gfc-v1";
    private static final String EMPTY = "-";

    private static final class StepResult {
        String output = "";
        String stdout = "";
        String stderr = "";
        String error = "";
    }

    private static final class Result {
        String status = "runner_failed";
        String failureStep = "runner";
        String intermediateWurcs = "";
        String iupacCondensed = "";
        String warning = "";
        String error = "";
    }

    private GfcBatchRunner() {}

    public static void main(String[] args) throws Exception {
        if (args.length != 1) {
            throw new IllegalArgumentException(
                "Usage: GfcBatchRunner <append-only-output.tsv>"
            );
        }
        File outputFile = new File(args[0]);
        File parent = outputFile.getAbsoluteFile().getParentFile();
        if (parent != null && !parent.exists() && !parent.mkdirs()) {
            throw new IllegalStateException("Unable to create " + parent);
        }
        try (
            BufferedReader input = new BufferedReader(
                new InputStreamReader(System.in, StandardCharsets.UTF_8)
            );
            PrintWriter output = new PrintWriter(
                new OutputStreamWriter(
                    new FileOutputStream(outputFile, true),
                    StandardCharsets.UTF_8
                )
            )
        ) {
            String line;
            int completed = 0;
            while ((line = input.readLine()) != null) {
                if (line.isEmpty()) continue;
                String[] fields = line.split("\\t", -1);
                if (fields.length != 4 || !PROTOCOL.equals(fields[0])) {
                    throw new IllegalArgumentException("Malformed request");
                }
                long requestId = Long.parseLong(fields[1]);
                Result result = convert(decode(fields[2]), decode(fields[3]));
                output.println(encodeResult(requestId, result));
                output.flush();
                completed++;
                if (completed % 1000 == 0) {
                    System.err.println(
                        "GlycanFormatConverter: " + completed + " rows complete"
                    );
                }
            }
        }
    }

    private static Result convert(String inputFormat, String sequence) {
        Result result = new Result();
        try {
            if ("GlycoCT".equals(inputFormat)) {
                StepResult first = invoke(
                    sequence.replace(" ", "\n"),
                    InputFormat.GLYCOCT,
                    OutputFormat.WURCS
                );
                result.intermediateWurcs = first.output;
                result.warning = combine(first.stdout, first.stderr);
                if (!first.error.isEmpty() || first.output.isEmpty()) {
                    result.status = "converter_failed";
                    result.failureStep = "glycoct_to_wurcs";
                    result.error = first.error.isEmpty()
                        ? swallowedFailure()
                        : first.error;
                    return result;
                }
                StepResult second = invoke(
                    first.output,
                    InputFormat.WURCS,
                    OutputFormat.IUPAC_CONDENSED
                );
                result.iupacCondensed = second.output;
                result.warning = combine(
                    result.warning,
                    combine(second.stdout, second.stderr)
                );
                if (!second.error.isEmpty() || second.output.isEmpty()) {
                    result.status = "converter_failed";
                    result.failureStep = "wurcs_to_iupac";
                    result.error = second.error.isEmpty()
                        ? swallowedFailure()
                        : second.error;
                    return result;
                }
            } else {
                InputFormat directFormat;
                if ("WURCS".equals(inputFormat)) {
                    directFormat = InputFormat.WURCS;
                } else if ("IUPAC-Extended".equals(inputFormat)) {
                    directFormat = InputFormat.IUPAC_EXTENDED;
                } else {
                    throw new IllegalArgumentException(
                        "Unsupported input format: " + inputFormat
                    );
                }
                StepResult direct = invoke(
                    sequence,
                    directFormat,
                    OutputFormat.IUPAC_CONDENSED
                );
                result.iupacCondensed = direct.output;
                result.warning = combine(direct.stdout, direct.stderr);
                if (!direct.error.isEmpty() || direct.output.isEmpty()) {
                    result.status = "converter_failed";
                    result.failureStep = "source_to_iupac";
                    result.error = direct.error.isEmpty()
                        ? swallowedFailure()
                        : direct.error;
                    return result;
                }
            }
            result.status = "converted";
            result.failureStep = "";
            return result;
        } catch (Throwable error) {
            result.error = stackTrace(error);
            return result;
        }
    }

    private static StepResult invoke(
        String sequence,
        InputFormat inputFormat,
        OutputFormat outputFormat
    ) {
        StepResult result = new StepResult();
        PrintStream originalOut = System.out;
        PrintStream originalErr = System.err;
        ByteArrayOutputStream stdout = new ByteArrayOutputStream();
        ByteArrayOutputStream stderr = new ByteArrayOutputStream();
        try (
            PrintStream capturedOut = new PrintStream(
                stdout,
                true,
                StandardCharsets.UTF_8.name()
            );
            PrintStream capturedErr = new PrintStream(
                stderr,
                true,
                StandardCharsets.UTF_8.name()
            )
        ) {
            System.setOut(capturedOut);
            System.setErr(capturedErr);
            ConverterPortal portal = new ConverterPortal();
            portal.inputSequence(sequence)
                .inputFormat(inputFormat)
                .outputFormat(outputFormat)
                .start();
            String converted = portal.getOutputSequence();
            result.output = converted == null ? "" : converted;
        } catch (Throwable error) {
            result.error = stackTrace(error);
        } finally {
            System.setOut(originalOut);
            System.setErr(originalErr);
            result.stdout = new String(
                stdout.toByteArray(),
                StandardCharsets.UTF_8
            );
            result.stderr = new String(
                stderr.toByteArray(),
                StandardCharsets.UTF_8
            );
        }
        return result;
    }

    private static String encodeResult(long requestId, Result result) {
        return String.join(
            "\t",
            PROTOCOL,
            Long.toString(requestId),
            result.status,
            result.failureStep,
            encode(result.intermediateWurcs),
            encode(result.iupacCondensed),
            encode(result.warning),
            encode(result.error)
        );
    }

    private static String combine(String left, String right) {
        if (left == null || left.isEmpty()) return right == null ? "" : right;
        if (right == null || right.isEmpty()) return left;
        return left + "\n" + right;
    }

    private static String swallowedFailure() {
        return "ConverterPortal returned an empty result; exception details " +
            "are unavailable because ConverterPortal.start() swallows them.";
    }

    private static String stackTrace(Throwable error) {
        java.io.StringWriter buffer = new java.io.StringWriter();
        error.printStackTrace(new PrintWriter(buffer));
        return buffer.toString();
    }

    private static String encode(String value) {
        if (value == null || value.isEmpty()) return EMPTY;
        return Base64.getEncoder().encodeToString(
            value.getBytes(StandardCharsets.UTF_8)
        );
    }

    private static String decode(String value) {
        if (EMPTY.equals(value)) return "";
        return new String(
            Base64.getDecoder().decode(value),
            StandardCharsets.UTF_8
        );
    }
}
