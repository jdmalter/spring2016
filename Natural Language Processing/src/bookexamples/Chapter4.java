package bookexamples;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import opennlp.tools.namefind.NameFinderME;
import opennlp.tools.namefind.NameSample;
import opennlp.tools.namefind.NameSampleDataStream;
import opennlp.tools.namefind.TokenNameFinderModel;
import opennlp.tools.tokenize.TokenizerME;
import opennlp.tools.tokenize.TokenizerModel;
import opennlp.tools.util.ObjectStream;
import opennlp.tools.util.PlainTextByLineStream;
import opennlp.tools.util.Span;

import com.aliasi.chunk.Chunk;
import com.aliasi.chunk.Chunker;
import com.aliasi.chunk.Chunking;
import com.aliasi.chunk.RegExChunker;
import com.aliasi.dict.DictionaryEntry;
import com.aliasi.dict.ExactDictionaryChunker;
import com.aliasi.dict.MapDictionary;
import com.aliasi.tokenizer.IndoEuropeanTokenizerFactory;
import com.aliasi.util.AbstractExternalizable;

/**
 * Chapter 4 examples and tests.
 * 
 * @author Jacob Malter learning from Natural Language Processing with Java
 *
 */
public class Chapter4 {

	private static String regularExpressionText //
	= "He left his email address (rgb@colorworks.com) and his "
			+ "phone number,800-555-1234. We believe his current address "
			+ "is 100 Washington Place, Seattle, CO 12345-1234. I "
			+ "understand you can also call at 123-555-1234 between "
			+ "8:00 AM and 4:30 most days. His URL is http://example.com "
			+ "and he was born on February 25, 1954 or 2/25/1954.";

	public static void javaRegex() {
		String phoneNumberRE = "\\d{3}-\\d{3}-\\d{4}";
		String urlRE = "\\b(https?|ftp|file|ldap)://[-A-Za-z0-9+&@#/%?=~_|!:,.;]*[-AZa-z0-9+&@#/%=~_|]";
		String zipRE = "[0-9]{5}(\\-?[0-9]{4})?";
		String emailRE = "[a-zA-Z0-9'._%+-]+@(?:[a-zA-Z0-9-]+\\.)+[a-zA-Z]{2,4}";
		String timeRE = "(([0-1]?[0-9])|([2][0-3])):([0-5]?[0-9])(:([0-5]?[0-9]))?";
		String dateRE = "((0?[13578]|10|12)(-|\\/)"
				+ "(([1-9])|(0[1-9])|([12])" + "([0-9]?)|(3[01]?))(-|\\/)"
				+ "((19)([2-9])(\\d{1})|(20)" + "([01])(\\d{1})|([8901])"
				+ "(\\d{1}))|(0?[2469]|11)(-|\\/)"
				+ "(([1-9])|(0[1-9])|([12])([0-" + "9]?)|(3[0]?))(-|\\/)((19)"
				+ "([2-9])(\\d{1})|(20)([01])" + "(\\d{1})|([8901])(\\d{1})))";
		Pattern pattern = Pattern.compile(phoneNumberRE + "|" + urlRE + "|"
				+ zipRE + "|" + emailRE + "|" + timeRE + "|" + dateRE);
		Matcher matcher = pattern.matcher(regularExpressionText);
		while (matcher.find())
			System.out.println(matcher.group() + " [" + matcher.start() + ":"
					+ matcher.end() + "]");
	}

	public static void lingpipeRegex() {
		Chunker chunker = new TimeRegexChunker();
		Chunking chunking = chunker.chunk(regularExpressionText);
		Set<Chunk> set = chunking.chunkSet();
		for (Chunk chunk : set)
			System.out.println("Type: "
					+ chunk.type()
					+ " Entity: ["
					+ regularExpressionText.substring(chunk.start(),
							chunk.end()) + "] Score: " + chunk.score());
	}

	@SuppressWarnings("serial")
	public static class TimeRegexChunker extends RegExChunker {
		private static final String TIME_RE = "(([0-1]?[0-9])|([2][0-3])):([0-5]?[0-9])(:([0-5]?[0-9]))?";
		private static final String CHUNK_TYPE = "time";
		private static final double CHUNK_SCORE = 1.0d;

		public TimeRegexChunker() {
			super(TIME_RE, CHUNK_TYPE, CHUNK_SCORE);
		}
	}

	private static String sentences[] = {
			"Joe was the last person to see Fred. ",
			"He saw him in Boston at McKenzie's pub at 3:00 where he "
					+ " paid $2.45 for an ale. ",
			"Joe wanted to go to Vermont for the day to visit a cousin who "
					+ "works at IBM, but Sally and he had to look for Fred" };

	public static void opennlpNER() {
		try {
			InputStream tokenStream = new FileInputStream(new File(
					"add file directory", "en-token.bin"));
			TokenizerModel tokenModel = new TokenizerModel(tokenStream);
			TokenizerME tokenizer = new TokenizerME(tokenModel);
			String modelNames[] = { "en-ner-person.bin", "en-ner-location.bin",
					"en-ner-organization.bin" };
			ArrayList<String> list = new ArrayList<>();
			for (String name : modelNames) {
				TokenNameFinderModel entityModel = new TokenNameFinderModel(
						new FileInputStream(
								new File("add file directory", name)));
				NameFinderME nameFinder = new NameFinderME(entityModel);
				for (int index = 0; index < sentences.length; index++) {
					String tokens[] = tokenizer.tokenize(sentences[index]);
					Span nameSpans[] = nameFinder.find(tokens);
					for (Span span : nameSpans)
						list.add("Sentence: " + index + " Span: "
								+ span.toString() + " Entity: "
								+ tokens[span.getStart()]);
				}
				for (String element : list)
					System.out.println(element);
			}
		} catch (Exception ex) {
			System.out.println("Exception");
		}
	}

	public static void lingpipeNER() {
		try {
			File modelFile = new File("add file directory",
					"ne-en-news-muc6.AbstractCharLmRescoringChunker");
			Chunker chunker = (Chunker) AbstractExternalizable
					.readObject(modelFile);
			for (int i = 0; i < sentences.length; ++i) {
				Chunking chunking = chunker.chunk(sentences[i]);
				System.out.println("Chunking=" + chunking);
			}
		} catch (IOException | ClassNotFoundException ex) {
			System.out.println("Exception");
		}
	}

	private static MapDictionary<String> dictionary;

	public static void lingpipeDictionary() {
		initializeDictionary();
		ExactDictionaryChunker dictionaryChunker = new ExactDictionaryChunker(
				dictionary, IndoEuropeanTokenizerFactory.INSTANCE, true, false);
		for (String sentence : sentences) {
			System.out.println("\nTEXT=" + sentence);
			Chunking chunking = dictionaryChunker.chunk(sentence);
			Set<Chunk> set = chunking.chunkSet();
			for (Chunk chunk : set) {
				System.out.println("Type: " + chunk.type() + " Entity: ["
						+ sentence.substring(chunk.start(), chunk.end())
						+ "] Score: " + chunk.score());
			}
		}
	}

	private static void initializeDictionary() {
		dictionary = new MapDictionary<String>();
		dictionary.addEntry(new DictionaryEntry<String>("Joe", "PERSON", 1.0));
		dictionary.addEntry(new DictionaryEntry<String>("Fred", "PERSON", 1.0));
		dictionary
				.addEntry(new DictionaryEntry<String>("Boston", "PLACE", 1.0));
		dictionary.addEntry(new DictionaryEntry<String>("pub", "PLACE", 1.0));
		dictionary
				.addEntry(new DictionaryEntry<String>("Vermont", "PLACE", 1.0));
		dictionary.addEntry(new DictionaryEntry<String>("IBM", "ORGANIZATION",
				1.0));
		dictionary
				.addEntry(new DictionaryEntry<String>("Sally", "PERSON", 1.0));
	}

	public static void opennlpTraining() {
		try (OutputStream modelOutputStream = new BufferedOutputStream(
				new FileOutputStream(new File("modelFile")))) {
			@SuppressWarnings("deprecation")
			ObjectStream<String> lineStream = new PlainTextByLineStream(
					new FileInputStream("en-ner-person.train"), "UTF-8");
			ObjectStream<NameSample> sampleStream = new NameSampleDataStream(
					lineStream);
			@SuppressWarnings("deprecation")
			TokenNameFinderModel model = NameFinderME.train("en", "person",
					sampleStream, Collections.<String, Object> emptyMap());
			model.serialize(modelOutputStream);
		} catch (IOException ex) {
			System.out.println("IO");
		}
	}

	/**
	 * Runs chapter 4 examples.
	 * 
	 * @param args
	 *            command line arguments not used in this program
	 */
	public static void main(String[] args) {
		javaRegex();
		System.out.println();
		lingpipeRegex();
		System.out.println();
		opennlpNER();
		System.out.println();
		lingpipeNER();
		System.out.println();
		lingpipeDictionary();
		System.out.println();
		opennlpTraining();
	}

}