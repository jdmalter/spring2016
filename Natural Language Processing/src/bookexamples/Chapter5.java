package bookexamples;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.ObjectInputStream;
import java.io.OutputStream;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;

import com.aliasi.classify.ConditionalClassification;
import com.aliasi.hmm.HiddenMarkovModel;
import com.aliasi.hmm.HmmDecoder;
import com.aliasi.tag.ScoredTagging;
import com.aliasi.tag.TagLattice;
import com.aliasi.tag.Tagging;
import com.aliasi.tokenizer.IndoEuropeanTokenizerFactory;
import com.aliasi.tokenizer.Tokenizer;
import com.aliasi.tokenizer.TokenizerFactory;

import opennlp.tools.chunker.ChunkerME;
import opennlp.tools.chunker.ChunkerModel;
import opennlp.tools.postag.MutableTagDictionary;
import opennlp.tools.postag.POSDictionary;
import opennlp.tools.postag.POSModel;
import opennlp.tools.postag.POSSample;
import opennlp.tools.postag.POSTaggerFactory;
import opennlp.tools.postag.POSTaggerME;
import opennlp.tools.postag.WordTagSampleStream;
import opennlp.tools.util.ObjectStream;
import opennlp.tools.util.PlainTextByLineStream;
import opennlp.tools.util.Sequence;
import opennlp.tools.util.Span;
import opennlp.tools.util.TrainingParameters;

/**
 * Chapter 5 examples and tests.
 * 
 * @author Jacob Malter learning from Natural Language Processing with Java
 *
 */
public class Chapter5 {

	private static String sentence[] = { "The", "voyage", "of", "the",
			"Abraham", "Lincoln", "was", "for", "a", "long", "time", "marked",
			"by", "no", "special", "incident." };

	private static String theSentence = "The voyage of the Abraham Lincoln was for a "
			+ "long time marked by no special incident.";

	public static void opennlpPOS() {
		try (InputStream modelIn = new FileInputStream(new File(
				"add file directory", "en-pos-maxent.bin"));) {
			POSModel model = new POSModel(modelIn);
			POSTaggerME tagger = new POSTaggerME(model);
			String tags[] = tagger.tag(sentence);
			for (int i = 0; i < sentence.length; i++)
				System.out.print(sentence[i] + "/" + tags[i] + " ");
			System.out.println();
			Sequence topSequences[] = tagger.topKSequences(sentence);
			for (int i = 0; i < topSequences.length; i++) {
				List<String> outcomes = topSequences[i].getOutcomes();
				double probabilities[] = topSequences[i].getProbs();
				for (int j = 0; j < outcomes.size(); j++) {
					System.out.printf("%s/%5.3f ", outcomes.get(j),
							probabilities[j]);
				}
				System.out.println();
			}
		} catch (IOException ex) {
			System.out.println("IO");
		}
	}

	public static void opennlpChunking() {
		try (InputStream posModelStream = new FileInputStream(
				"add file directory\\en-pos-maxent.bin");
				InputStream chunkerStream = new FileInputStream(
						"add file directory\\en-chunker.bin");) {
			POSModel model = new POSModel(posModelStream);
			POSTaggerME tagger = new POSTaggerME(model);
			String tags[] = tagger.tag(sentence);
			for (int i = 0; i < tags.length; i++)
				System.out.print(sentence[i] + "/" + tags[i] + " ");
			System.out.println();
			ChunkerModel chunkerModel = new ChunkerModel(chunkerStream);
			ChunkerME chunkerME = new ChunkerME(chunkerModel);
			String result[] = chunkerME.chunk(sentence, tags);
			for (int i = 0; i < result.length; i++)
				System.out.println("[" + sentence[i] + "]" + result[i]);
			Span spans[] = chunkerME.chunkAsSpans(sentence, tags);
			for (Span span : spans) {
				System.out.print("Type: " + span.getType() + " - " + " Begin: "
						+ span.getStart() + " End: " + span.getEnd()
						+ " Length: " + span.length() + "  [ ");
				for (int j = span.getStart(); j < span.getEnd(); j++)
					System.out.print(sentence[j] + " ");
				System.out.println("]");
			}
		} catch (IOException ex) {
			System.out.println("IO");
		}
	}

	public static void opennlpDictionary() {
		try (InputStream modelIn = new FileInputStream(new File(
				"add file directory", "en-pos-maxent.bin"));) {
			POSModel model = new POSModel(modelIn);
			POSTaggerFactory posTaggerFactory = model.getFactory();
			MutableTagDictionary tagDictionary = (MutableTagDictionary) posTaggerFactory
					.getTagDictionary();
			String tags[] = tagDictionary.getTags("force");
			for (String tag : tags)
				System.out.print("/" + tag);
			System.out.println();
			String oldTags[] = tagDictionary.put("force", "newTag");
			for (String tag : oldTags)
				System.out.print("/" + tag);
			System.out.println();
			String newTags[] = new String[tags.length + 1];
			for (int i = 0; i < tags.length; i++)
				newTags[i] = tags[i];
			newTags[tags.length] = "newTag";
			oldTags = tagDictionary.put("force", newTags);
			tags = tagDictionary.getTags("force");
			for (String tag : tags)
				System.out.print("/" + tag);
			System.out.println();
		} catch (IOException ex) {
			System.out.println("IO");
		}
	}

	public static void opennlpCreateDictionary() {
		try (InputStream dictionaryIn = new FileInputStream(new File(
				"dictionary.xml"));) {
			POSDictionary dictionary = POSDictionary.create(dictionaryIn);
			Iterator<String> iterator = dictionary.iterator();
			while (iterator.hasNext()) {
				String entry = iterator.next();
				String tags[] = dictionary.getTags(entry);
				System.out.print(entry + " ");
				for (String tag : tags)
					System.out.print("/" + tag);
				System.out.println();
			}
		} catch (IOException ex) {
			System.out.println("IO");
		}
	}

	public static void lingpipeTagger() {
		try (FileInputStream inputStream = new FileInputStream(
				"add file directory\\pos-en-general-brown.HiddenMarkovModel");
				ObjectInputStream objectStream = new ObjectInputStream(
						inputStream);) {
			HiddenMarkovModel hmm = (HiddenMarkovModel) objectStream
					.readObject();
			HmmDecoder decoder = new HmmDecoder(hmm);
			TokenizerFactory TOKENIZER_FACTORY = IndoEuropeanTokenizerFactory.INSTANCE;
			char charArray[] = theSentence.toCharArray();
			Tokenizer tokenizer = TOKENIZER_FACTORY.tokenizer(charArray, 0,
					charArray.length);
			String tokens[] = tokenizer.tokenize();
			List<String> tokenList = Arrays.asList(tokens);
			Tagging<String> tagString = decoder.tag(tokenList);
			for (int i = 0; i < tagString.size(); ++i)
				System.out.print(tagString.token(i) + "/" + tagString.tag(i)
						+ " ");
			System.out.println();
			String newSentence[] = { "Bill", "used", "the", "force", "to",
					"force", "the", "manager", "to", "tear", "the", "bill",
					"in", "two." };
			tokenList = Arrays.asList(newSentence);
			int maxResults = 5;
			Iterator<ScoredTagging<String>> iterator = decoder.tagNBest(
					tokenList, maxResults);
			while (iterator.hasNext()) {
				ScoredTagging<String> scoredTagging = iterator.next();
				System.out.printf("Score: %7.3f   Sequence: ",
						scoredTagging.score());
				for (int i = 0; i < tokenList.size(); i++)
					System.out.print(scoredTagging.token(i) + "/"
							+ scoredTagging.tag(i) + " ");
				System.out.println();
			}
			TagLattice<String> lattice = decoder.tagMarginal(tokenList);
			for (int index = 0; index < tokenList.size(); index++) {
				ConditionalClassification classification = lattice
						.tokenClassification(index);
				System.out.printf("%-8s", tokenList.get(index));
				for (int i = 0; i < 4; ++i) {
					double score = classification.score(i);
					String tag = classification.category(i);
					System.out.printf("%7.3f/%-3s ", score, tag);
				}
				System.out.println();
			}
		} catch (IOException ex) {
			System.out.println("IO");
		} catch (ClassNotFoundException ex) {
			System.out.println("ClassNotFound");
		}
	}

	@SuppressWarnings("deprecation")
	public static void opennlpTrainingModel() {
		POSModel model = null;
		try (InputStream dataIn = new FileInputStream("sample.train");) {
			ObjectStream<String> lineStream = new PlainTextByLineStream(dataIn,
					"UTF-8");
			ObjectStream<POSSample> sampleStream = new WordTagSampleStream(
					lineStream);
			model = POSTaggerME.train("en", sampleStream,
					TrainingParameters.defaultParams(), null, null);
			try (OutputStream modelOut = new BufferedOutputStream(
					new FileOutputStream(new File("en_pos_verne.bin")));) {
				model.serialize(modelOut);
			} catch (IOException ex) {
				System.out.println("IO in IO");
			}
		} catch (IOException ex) {
			System.out.println("IO");
		}
	}

	/**
	 * Runs chapter 5 examples.
	 * 
	 * @param args
	 *            command line arguments not used in this program
	 */
	public static void main(String[] args) {
		opennlpPOS();
		System.out.println();
		opennlpChunking();
		System.out.println();
		opennlpDictionary();
		System.out.println();
		opennlpCreateDictionary();
		System.out.println();
		lingpipeTagger();
		System.out.println();
		opennlpTrainingModel();
	}

}