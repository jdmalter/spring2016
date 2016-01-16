package bookexamples;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.io.Reader;
import java.io.StringReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;
import java.util.Properties;

import opennlp.tools.cmdline.postag.POSModelLoader;
import opennlp.tools.namefind.NameFinderME;
import opennlp.tools.namefind.TokenNameFinderModel;
import opennlp.tools.postag.POSModel;
import opennlp.tools.postag.POSTaggerME;
import opennlp.tools.tokenize.SimpleTokenizer;
import opennlp.tools.tokenize.Tokenizer;
import opennlp.tools.tokenize.TokenizerME;
import opennlp.tools.tokenize.TokenizerModel;
import opennlp.tools.tokenize.WhitespaceTokenizer;
import opennlp.tools.util.Span;

import com.aliasi.tokenizer.IndoEuropeanTokenizerFactory;

import edu.stanford.nlp.ling.CoreLabel;
import edu.stanford.nlp.ling.HasWord;
import edu.stanford.nlp.pipeline.Annotation;
import edu.stanford.nlp.pipeline.StanfordCoreNLP;
import edu.stanford.nlp.process.CoreLabelTokenFactory;
import edu.stanford.nlp.process.DocumentPreprocessor;
import edu.stanford.nlp.process.PTBTokenizer;

/**
 * Chapter 1 examples and tests.
 * 
 * @author Jacob Malter learning from Natural Language Processing with Java
 *
 */
public class Chapter1 {

	public static void apacheOpenNLP() {
		try (InputStream is = new FileInputStream(new File(
				"add file directory", "en-token.bin"))) {
			TokenizerModel model = new TokenizerModel(is);
			Tokenizer tokenizer = new TokenizerME(model);
			String[] tokens = tokenizer.tokenize("He lives at 1511 W. "
					+ "Randolph.");
			for (String a : tokens) {
				System.out.print("[" + a + "] ");
			}
			System.out.println();
		} catch (FileNotFoundException e) {
			System.out.println("FileNotFoundException thrown");
		} catch (IOException ex) {
			System.out.println("IOException thown");
		}
	}

	public static void standfordNLP() {
		PTBTokenizer<CoreLabel> ptb = new PTBTokenizer<>(new StringReader(
				"He lives at 1511 W. Randolph."), new CoreLabelTokenFactory(),
				null);
		while (ptb.hasNext()) {
			System.out.print("(" + ptb.next() + ") ");
		}
		System.out.println();
	}

	public static void lingPipe() {
		List<String> tokenList = new ArrayList<String>();
		List<String> whiteList = new ArrayList<String>();
		String text = "A sample sentence processed \nby \tthe "
				+ "LingPipe tokenizer";
		com.aliasi.tokenizer.Tokenizer tokenizer = IndoEuropeanTokenizerFactory.INSTANCE
				.tokenizer(text.toCharArray(), 0, text.length());
		tokenizer.tokenize(tokenList, whiteList);
		for (String element : tokenList) {
			System.out.print(element + " ");
		}
		System.out.println();
	}

	/**
	 * Related to Chapter 2, Finding Parts of Text.
	 */
	public static void javaDelimiters() {
		String text = "Mr. Smith went to 123 Washington avenue.";
		String[] tokens = text.split("\\s+");
		for (String token : tokens) {
			System.out.print("<" + token + "> ");
		}
		System.out.println();
	}

	/**
	 * Related to Chapter 3, Finding Sentences.
	 */
	public static void sentenceFinder() {
		String paragraph = "The first sentence. The second sentence.";
		Reader reader = new StringReader(paragraph);
		DocumentPreprocessor documentPreprocessor = new DocumentPreprocessor(
				reader);
		List<String> sentenceList = new LinkedList<String>();
		for (List<HasWord> element : documentPreprocessor) {
			StringBuilder sentence = new StringBuilder();
			List<HasWord> hasWordList = element;
			for (HasWord token : hasWordList) {
				sentence.append(token).append(" ");
			}
			sentenceList.add(sentence.toString());
		}
		for (String sentence : sentenceList) {
			System.out.println(sentence);
		}
	}

	/**
	 * Related to Chapter 4, Finding People and Things.
	 */
	public static void peopleFinder() {
		try {
			String[] sentences = { "Tim was a good neighbor. Perhaps not as good as a Bob "
					+ "Haywood, but still pretty good. Of course Mr. Adam "
					+ "took the cake!" };
			Tokenizer tokenizer = SimpleTokenizer.INSTANCE;
			TokenNameFinderModel model = new TokenNameFinderModel(new File(
					"add file directory", "en-ner-person.bin"));
			NameFinderME finder = new NameFinderME(model);
			for (String sentence : sentences) {
				String[] tokens = tokenizer.tokenize(sentence);
				Span[] nameSpans = finder.find(tokens);
				System.out.println(Arrays.toString(Span.spansToStrings(
						nameSpans, tokens)));
			}
		} catch (IOException e) {
			System.out.println();
		}
	}

	/**
	 * Related to Chapter 5, Detecting Parts of Speech.
	 */
	public static void partOfSpeechFinder() {
		POSModel model = new POSModelLoader().load(new File(
				"add file directory", "en-pos-maxent.bin"));
		POSTaggerME tagger = new POSTaggerME(model);
		String sentence = "POS processing is useful for enhancing the "
				+ "quality of data sent to other elements of a pipeline.";
		String[] tokens = WhitespaceTokenizer.INSTANCE.tokenize(sentence);
		String[] tags = tagger.tag(tokens);
		for (int i = 0; i < tokens.length; i++) {
			System.out.print(tokens[i] + "[" + tags[i] + "] ");
		}
	}

	// There are no examples for Chapter 6, Classifying Text and Documents.

	/**
	 * Related to Chapter 7, Using a Parser to Extract Relationships. However,
	 * StandfordCoreNLP class cannot instantiate.
	 */
	@SuppressWarnings("unused")
	private static void relationshipFinder() {
		Properties properties = new Properties();
		properties.put("annotators", "tokenize, ssplit, parse");
		StanfordCoreNLP pipeline = new StanfordCoreNLP(properties);
		Annotation annotation = new Annotation(
				"The meaning and purpose of life is plain to see.");
		pipeline.annotate(annotation);
		pipeline.prettyPrint(annotation, System.out);
	}

	// There are no examples for Chapter 8, Combined Approaches.

	/**
	 * Runs chapter 1 examples.
	 * 
	 * @param args
	 *            command line arguments not used in this program
	 */
	public static void main(String[] args) {
		apacheOpenNLP();
		standfordNLP();
		lingPipe();
		javaDelimiters();
		sentenceFinder();
		peopleFinder();
		partOfSpeechFinder();
	}

}