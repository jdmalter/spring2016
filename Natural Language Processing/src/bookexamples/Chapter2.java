package bookexamples;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.StringReader;
import java.text.BreakIterator;
import java.util.ArrayList;
import java.util.List;
import java.util.Scanner;

import edu.stanford.nlp.ling.CoreLabel;
import edu.stanford.nlp.process.CoreLabelTokenFactory;
import edu.stanford.nlp.process.PTBTokenizer;
import opennlp.tools.tokenize.SimpleTokenizer;
import opennlp.tools.tokenize.Tokenizer;
import opennlp.tools.tokenize.TokenizerME;
import opennlp.tools.tokenize.TokenizerModel;
import opennlp.tools.tokenize.WhitespaceTokenizer;

/**
 * Chapter 2 examples and tests.
 * 
 * @author Jacob Malter learning from Natural Language Processing with Java
 *
 */
public class Chapter2 {

	private static String paragraph = "Let's pause, \nand then " + "reflect.";

	public static void javaTokenizers() {
		Scanner scanner = new Scanner("Let's pause, and then " + "reflect.");
		scanner.useDelimiter("[ ,.]");
		List<String> list = new ArrayList<String>();
		while (scanner.hasNext()) {
			list.add(scanner.next());
		}
		for (String token : list) {
			System.out.print(token + "_");
		}
		System.out.println();
		scanner.close();
	}

	public static void javaStringSplit() {
		String text = "Mr. Smith went to 123 Washington avenue.";
		String[] tokens = text.split("\\s+");
		for (String token : tokens) {
			System.out.print("[" + token + "]");
		}
		System.out.println();
	}

	public static void javaBreakIterator() {
		BreakIterator wordIterator = BreakIterator.getWordInstance();
		String text = "Let's pause, and then reflect.";
		wordIterator.setText(text);
		int boundary = wordIterator.first();
		while (boundary != BreakIterator.DONE) {
			int begin = boundary;
			System.out.print(boundary + "-");
			boundary = wordIterator.next();
			int end = boundary;
			if (end == BreakIterator.DONE)
				break;
			System.out.print(boundary + " [" + text.substring(begin, end)
					+ "];");
		}
		System.out.println();
	}

	public static void openNLPTokenizer() {
		SimpleTokenizer simpleTokenizer = SimpleTokenizer.INSTANCE;
		String[] tokens = simpleTokenizer.tokenize(paragraph);
		for (String token : tokens) {
			System.out.print("[" + token + "] ");
		}
		System.out.println();
		WhitespaceTokenizer whiteSpaceTokenizer = WhitespaceTokenizer.INSTANCE;
		String[] whitespaceTokens = whiteSpaceTokenizer.tokenize(paragraph);
		for (String token : whitespaceTokens) {
			System.out.print("[" + token + "] ");
		}
		System.out.println();
	}

	public static void maximumEntropyTokenizer() {
		try {
			InputStream modelInputStream = new FileInputStream(new File(
					"add file directory", "en-token.bin"));
			TokenizerModel model = new TokenizerModel(modelInputStream);
			Tokenizer tokenizer = new TokenizerME(model);
			String[] tokens = tokenizer.tokenize(paragraph);
			for (String token : tokens) {
				System.out.print("(" + token + ") ");
			}
			System.out.println();
		} catch (IOException e) {
			System.out.println("IOException thrown.");
		}
	}

	public static void standfordNLP() {
		CoreLabelTokenFactory ctf = new CoreLabelTokenFactory();
		PTBTokenizer<CoreLabel> ptb = new PTBTokenizer<>(new StringReader(
				paragraph), ctf, "invertible=true");
		while (ptb.hasNext()) {
			CoreLabel cl = ptb.next();
			System.out.print(cl.originalText() + " [" + cl.beginPosition()
					+ "-" + cl.endPosition() + "];");
		}
		System.out.println();
	}

	/**
	 * Runs chapter 2 examples.
	 * 
	 * @param args
	 *            command line arguments not used in this program
	 */
	public static void main(String[] args) {
		javaTokenizers();
		javaStringSplit();
		javaBreakIterator();
		openNLPTokenizer();
		maximumEntropyTokenizer();
		standfordNLP();
	}

}